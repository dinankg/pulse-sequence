%% 3D GRE pulse sequence coded with PulSeq
%{
Dinank Gupta 
University of Michigan, June 2024
TODO: Change RF to SLR/SPSP
%}

clc
clear

%% Constants:
gamT = 4.2576e7;   % Hz/Tesla
gamG = gamT/1e4;    % Hz/Gauss

%% Flag to save ktraj
save_ktraj=1;
%%
%% Basic Scan Parameters

seq_params.scanner = 'inside'; %For use at Michigan
seq_params.trig = 0;
seq_params.TR = 40e-3;
%TE is set to minTE for now.
seq_params.ntp = 5;
seq_params.nDummyLoops = 2;
seq_params.nshot_spiral = 12;
seq_params = getparams(seq_params); % Important scan info here
seq_params.oversamp=100;%Fully sample these many spiral samples.
seq_params.scanner_type = 'philips'; % or 'ge' or 'seimens'
%Option to turn fatsat on/off
seq_params.dofatsat = 1;

sys=mr.opts('maxGrad',seq_params.maxgrad,'gradUnit', ...
    'mT/m','maxSlew',seq_params.maxslew, 'slewUnit', 'T/m/s');
sys.adcDeadTime=1e-5;
% sys.rfRasterTime=1e-5;
seq = mr.Sequence(sys);

%% Fatsat pulse: Copied from Pulseq demo
% Create fat-sat pulse
sat_freq=seq_params.fatChemShift*seq_params.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
rf_fs.deadTime=100e-6;% 100us deadtime (mostly needed for GE)
rf_fs.delay =rf_fs.deadTime;

gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/(seq_params.fov(1)/seq_params.N(1))*6); % 6 cycles of phase over x res.
%%
% RF pulse
seq_params.rf_slice_fact=0.8; % Excite a smaller slice than prescribed
[rf, gz_rf,gzReph] = mr.makeSincPulse(seq_params.alpha*pi/180, 'Duration', 3e-3, ...
    'SliceThickness', seq_params.slicethickness*seq_params.rf_slice_fact, 'apodization', 0.42, ...
    'timeBwProduct', 6, 'system', sys);
gzPreph = gzReph;
rf.delay = rf.deadTime;
gz_rf.delay = ceil((mr.calcDuration(gzReph))/sys.gradRasterTime)*sys.gradRasterTime;
rf.delay = gz_rf.delay+gz_rf.riseTime; %-5e-7
gzReph.delay =  mr.calcDuration(gz_rf);
gzcomb = mr.addGradients({gzPreph,gz_rf,gzReph},'system',sys);
%% Spiral segment
% create spiral readout with multiple shot
[gx1,gy1,t,spiral_readout_length] = makevdspiral(seq_params.fov(1)*100,seq_params.N(1), ...
    seq_params.nshot_spiral,seq_params.oversamp,seq_params.maxgrad/10, ...
    0.9*seq_params.maxslew ,sys.gradRasterTime);
%making it divisible by 10
gx1 = [gx1;zeros(10-rem(length(gx1),10),size(gx1,2))];
gy1 = [gy1;zeros(10-rem(length(gy1),10),size(gy1,2))];

gx1_Hzcm = gx1*gamG*100;gy1_Hzcm = gy1*gamG*100; %Converting to Hz/m
% Making struct with all spiral shots.
gx = cell(1,seq_params.nshot_spiral);
gy = cell(1,seq_params.nshot_spiral);
% last_sample_time = sys.gradRasterTime*(size(gx1,1)-1);
for i = 1:seq_params.nshot_spiral
    gx{i} = mr.makeArbitraryGrad('x',gx1_Hzcm(:,i),sys,'delay',.1e-3);
    gy{i}  = mr.makeArbitraryGrad('y',gy1_Hzcm(:,i),sys,'delay',.1e-3);
end

%% kz-encoding segment
dkz = 1/seq_params.fov(3); %1/m
res = seq_params.fov(3)/seq_params.N(3);%m
seq_params.kzmax = 1/res;%1/m
z_enc_trap = mr.makeTrapezoid('z','Area',seq_params.kzmax,'system',sys);

% Grads will go from -Gz_area/2 to Gz_area/2
zp = (-(seq_params.N(3)-1)/2:(seq_params.N(3)-1)/2)/((seq_params.N(3)-1)/2)/2; %zpe scaling for each TR
%List of gradient area:
seq_params.zp_scale = zp;%/zp(1); %Scaling of the kz traps

%% ADC
adc = mr.makeAdc(spiral_readout_length,sys,'Dwell',sys.gradRasterTime,'delay',.1e-3);
%% Spoiler
g_spoilx = mr.makeTrapezoid('x', ...
    'Area', 1/(seq_params.fov(1)/seq_params.N(1))*2, ...  % spoil with 2x cycles per voxel /10 to convert to m
    'system',sys);
g_spoilz= mr.makeTrapezoid('z', ...
    'Area', 1/(seq_params.fov(3)/seq_params.N(3))*2, ...  % spoil with 2x cycles per voxel
    'system',sys);

%% Define delays
TRmin = mr.calcDuration(gzcomb)+mr.calcDuration(gx{1}) + 2*mr.calcDuration(z_enc_trap)+mr.calcDuration(g_spoilz);
TRmin = TRmin + seq_params.dofatsat*mr.calcDuration(gz_fs);
delayTR = ceil((seq_params.TR - TRmin)/seq.gradRasterTime)*seq.gradRasterTime;

%% Putting it all together

ndisdaq = seq_params.nDummyLoops;
adc_flag=0;
for itp = 0:seq_params.ntp % 0th timepoint will serve as disdaq.
    for nd=1:ndisdaq % n_disdaq will be set to 1 after adc is on.
        if(adc_flag==0) %No adc
            segmentID = 0;
        elseif(adc_flag==1)
            segmentID=seq_params.nshot_spiral;
        end
        for ishots = 1:seq_params.nshot_spiral

            for ikz = 1:seq_params.N(3)

                segmentID = ishots+1 + (adc_flag*seq_params.nshot_spiral);
                % transmit and receive phase cycling
                rf.phaseOffset = rf.phaseOffset+deg2rad(117);
                adc.phaseOffset = adc.phaseOffset+deg2rad(117);

                if(seq_params.dofatsat)
                    %fat-sat and slab select
                    seq.addBlock(rf_fs,gz_fs, mr.makeLabel('SET', 'TRID', segmentID))%
                    seq.addBlock(rf,gzcomb)
                else
                    % Only slab select
                    seq.addBlock(rf,gzcomb, mr.makeLabel('SET', 'TRID', segmentID))
                end

                %kz-encode
                seq.addBlock(mr.scaleGrad(z_enc_trap,seq_params.zp_scale(ikz)));
                if adc_flag
                    seq.addBlock(gx{ishots},gy{ishots},adc)
                else
                    seq.addBlock(gx{ishots},gy{ishots})
                end
                %kz-encode balance
                seq.addBlock(mr.scaleGrad(z_enc_trap,-1*seq_params.zp_scale(ikz)));
                % spoiler
                seq.addBlock(g_spoilx,g_spoilz);

                seq.addBlock(mr.makeDelay(delayTR));

            end
        end
        if(nd==ndisdaq && adc_flag==0)% After the disdaq loop ends once, set adc to 1.
            adc_flag=1;
            ndisdaq=1;%After disdaqs are done, set the ndisdaq to 1
            break
        end
    end

end

%%
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.plot('timerange',[0 seq_params.TR]);
%% Output for execution
seq.setDefinition('FOV', seq_params.fov);
seq.setDefinition('Name', '3DGRE');
seq.write('3dgre.seq');

%% Seeing in GE format to make sure units look ok
sysGE = toppe.systemspecs('maxGrad', 10, ... % G/cm
    'maxSlew', 20, ... % G/cm/ms
    'maxRF', 0.2, ... % Gauss. Must be >= peak RF in sequence
    'rfDeadTime',100,...
    'adcDeadTime', 0, ... % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ... % RF/gradient delay (us)
    'psd_grd_wait', 156); % ADC/gradient delay (us)

seq2ge('3dgre.seq',sysGE,'3dgre.tar');

system('tar -xvf 3dgre.tar');

figure; toppe.plotseq(sysGE,'timeRange',[0 seq_params.TR]); % Plot 1 frame

%% Plotting multiple TRs:
%{
for k=1:5:300
    toppe.plotseq(sysGE,'timeRange',[k*seq_params.TR, (k+1)*seq_params.TR]);
    drawnow
end
%}
if(save_ktraj)
    dt = adc.dwell;%sec
    ktraj = cat(3,cumtrapz(gx1),cumtrapz(gy1))*dt*gamG*100; %1/m % klen x nshots x 2
    ktraj = permute(ktraj,[1,3,2]);% klen x 2 x nshots
    n_adc = floor(adc.numSamples/4)*4;
    ktraj=ktraj(1:n_adc,:,:); % Cropping only the acquired part
    ks = ktraj.*(2*pi*seq_params.fov(1)/seq_params.N(1));
    %% Sorting the kdata and ktrah according to shots
    ks = squeeze(ks);kspace=zeros(size(ks,1),3,seq_params.N(3),seq_params.nshot_spiral);
    n=0;
    for nt = 1: seq_params.ntp
        for ns = 1:seq_params.nshot_spiral
            for kz = 1:seq_params.N(3)
                kz_tmp = seq_params.zp_scale(kz)*repmat(seq_params.kzmax ,[length(ks),1]);
                ks(:,3,ns) = kz_tmp*(2*pi*seq_params.fov(3)/seq_params.N(3));
                kspace(:,:,kz,ns) = ks(:,:,ns);
                n=n+1;
            end
        end
    end
    % kspace size is klen x 3 x kz x nshots
    % Assumes same ktraj across timepoints
    sp_number=0
    while true
        sp_number = sp_number + 1
        sp_name = ['spiral_ktraj',num2str(sp_number),'.mat'];

        if ~isfile(sp_name)
            sprintf('Saving kspace traj as spiral_ktraj%1d...',sp_number)
            save(sp_name,'kspace','seq_params')
            break
        end
    end
end
delete('*.mod')

%%
% if ACTUAL_SCAN
%     sprintf('This is used internally at fMRI lab at UMICH.')
%     exp_number = 0;
%     while true
%         exp_number = exp_number + 1;
%         exp_name = ['experiments/',datestr(now,'yyyy-mm-dd_'),seq_params.scanner,'_exp',num2str(exp_number)];
%         if ~isfolder(exp_name)
%             mkdir(exp_name)
%             break
%         end
%     end
%     [status] = copyfile('3dgre.tar', ['./',exp_name,'/3dgre.tar']);
%     [status] = copyfile('3dgre.seq', ['./',exp_name,'/3dgre.seq']);
%     system(append('tar -xvf ','./',exp_name,'/','3dgre.tar',...
%         ' -C ',exp_name,'/'));
%
%
%     save('variables')
%     [status] = copyfile('variables.mat', ['./',exp_name,'/variables.mat']);
%
%
%     toppe.utils.coppe('target',seq_params.scanner,'cv',788,'use_pw',1);
% end

