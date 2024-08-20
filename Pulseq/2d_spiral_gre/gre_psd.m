% Make 2D GRE with spiral acquisition using Pulseq.

% Also needs rf_tools and MIRT packages.
% Can also get Larson's Spectral-Spatial-RF-Pulse-Design package to estimate slice profile
% Dinank Gupta Aug 2024
% University of Michigan, Ann Arbor

clc
clear
%% Constants
gamT = 4.2576e7;   % Hz/Tesla
gamG = gamT/1e4;    % Hz/Gauss
seq.gmax_Gcm = 5; %G/cm
seq.smax_mTms = 200;%mT/m/ms
save_ktraj = 1;
%% Scanner Parameters:
scanner = 'GE';
if strcmp(scanner,'Siemens')
    sys = mr.opts('MaxGrad', seq.gmax_Gcm *10, 'GradUnit', 'mT/m', ...
        'MaxSlew', seq.smax_Gcmms, 'SlewUnit', 'T/m/s', ...
        'B0', 123/128*3, ...
        'rfDeadTime', 100e-6, ...
        'rfRingdownTime', 10e-6, ...
        'adcDeadTime', 10e-6);
elseif strcmp(scanner,'GE')
    sys = mr.opts('MaxGrad', seq.gmax_Gcm *10, 'GradUnit', 'mT/m', ...
        'MaxSlew', seq.smax_mTms, 'SlewUnit', 'T/m/s', 'B0', 3.0, ...
        'rfRasterTime',4e-6,'rfDeadTime', 100e-6,'rfRingdownTime', 60e-6, ...
        'gradRasterTime',4e-6,'blockDurationRaster',4e-6,...
        'adcRasterTime',4e-6,'adcDeadTime', 0e-6);
end
psd = mr.Sequence(sys); %#ok<*DPSD>  % Create a new sequence object
seq.rf_dt = psd.rfRasterTime;

%% Sequence Parameters:
% Field of view, Resolution
seq.FOV = 24*1e-2; % Field of view (m) in x-y
seq.nshot = 20; % # of spiral rotations
seq.npix = 256; % # of voxels in x-y
seq.slthick = 0.003;     % Slice thickness in m
seq.nsli=15; % set it to number of slices you want in fully-sampled case.
seq.slspacing = 0.001; % Spacing between slices in m

seq.TE = 5;%ms
seq.TR = 1000;%ms
seq.TR_eff = seq.TR/seq.nsli; % Effective scan TR (for interleaving)

seq.ntp = 2; %time points for dynamic scan
seq.spiraldir = 1;     % 1 is spiral out, 2 is spiral in

% Using Ernst angle for T1 of 1000ms
seq.fa=1*acosd(exp(-seq.TR/1000));
%% Flags:
seq.addfatsat = 1; % To have fatsat or not.
%% Fatsat pulse:
seq.fat_dur = 2 * 1e-3;%s
seq.fat_tbw = 3;%TBW
seq.fat_offset = -440;%Hz

rf_fat = dzrf(round(seq.fat_dur/seq.rf_dt),seq.fat_tbw,'ex','ls',0.001,0.001); % pm excitation pulse
rf_fat = rf_fat/sum(rf_fat);
rf_fatsat = mr.makeArbitraryRf(rf_fat, pi/2, 'system', sys);
g_fatsat= mr.makeTrapezoid('z','Area', 1/(seq.FOV/seq.npix)*2, ...  % spoil with 2x cycles per voxel
    'Delay',rf_fatsat.deadTime+rf_fatsat.shape_dur+rf_fatsat.ringdownTime,...
    'system',sys);

%% Excitation pulse:
seq.ex_dur = 3 *1e-3;      % RF pulse duration in ms
seq.ex_tbw = 6;           % Time bandwidth for SLR pulse
if seq.fa < 20, seq.ex_type = 'st'; else, seq.ex_type = 'ex'; end
seq.ex_ftype = 'pm';

[rf,grad] =  slr_pulse(seq.slthick,seq.ex_tbw,seq.ex_dur,seq.ex_type,seq.ex_ftype,sys);
rf_ex = mr.makeArbitraryRf([rf], deg2rad(seq.fa), 'system', sys);
grad_ex = mr.makeArbitraryGrad('z',[grad*gamG*100],sys);
[f,z,m] = ss_plot(grad_ex.waveform/(gamG*100), rf_ex.signal/1e4, 4e-6,[],4,1000,[],[0] ,[]);
vline(-seq.slthick/2*100),vline(seq.slthick/2*100)

% Setting up list of frequency offsets
seq.rf_bw = (seq.ex_tbw/seq.ex_dur)/seq.slthick; %BW of RF pulses in Hz/m
seq.sl_locs = linspace(-floor(seq.nsli/2),floor(seq.nsli/2),seq.nsli)*(seq.slthick+seq.slspacing);% Center of slice locs m.
seq.sl_freqs = seq.rf_bw*seq.sl_locs; % Frequency for slice sel in Hz.
%% Readout
seq.oversamp=200; % Number of samples to acquire at fully sampled rate.
[gx,gy,t,npts] = makevdspiral2(seq.FOV*100,seq.npix,seq.nshot,seq.oversamp,seq.gmax_Gcm,0.7*seq.smax_mTms);
gx = gx*gamG*100;gy = gy*gamG*100; %Converting to Hz/m
% gz=0*gx;
% paramsint16(2) = npts;% p
seq.spiral_pts = npts;

for k=1:size(gx,2)
    gx_read{k} = mr.makeArbitraryGrad('x',gx(:,k),sys);
    gy_read{k} = mr.makeArbitraryGrad('y',gy(:,k),sys);
end
% List of rotations for spirals:
seq.shotrot = linspace(0,2*pi-(2*pi)/seq.nshot,seq.nshot);
seq.rot_list = cat(2,repmat([0 0 1],seq.nshot,1),seq.shotrot.');
seq.rotmat = axang2rotm(seq.rot_list);
seq.readout_grad_mT_m = cat(3,gx/(gamG*10),gy/(gamG*10)); %mT/m
%% ADC CHECK FOR CORRECTNESS!!!
num_adc_samp = npts;ceil(gx_read{1}.shape_dur/sys.adcRasterTime);
%Add delay associated to the shootout trapezoid asociated with spiral in.
adc=mr.makeAdc(num_adc_samp, 'Dwell', sys.adcRasterTime, 'Delay', sys.adcDeadTime,'system',sys);
%% Spoiler
g_spoilx = mr.makeTrapezoid('x', ...
    'Area', 1/(seq.FOV/seq.npix)*6, ...  % spoil with 6x cycles per voxel /10 to convert to m
    'system',sys);
g_spoilz= mr.makeTrapezoid('z', ...
    'Area', 1/(seq.FOV/seq.npix)*6, ...  % spoil with 6x cycles per voxel
    'system',sys);

%% Getting the timing right:
t_current = seq.addfatsat*(rf_fatsat.shape_dur + mr.calcDuration(g_fatsat))+...
    rf_ex.shape_dur +...
    gx_read{1}.shape_dur+...
    mr.calcDuration(g_spoilz);
delayTR = round(seq.TR_eff*1e-3 - t_current,3);
%% Combining everything:

n_disdaq = 15;% TRs after which daq will turn on.
adc_flag = 0;
sp_shot=1;% flag for spiral
%
% Making scan loop
for itp = 0:seq.ntp
    % itp
    for ishot = 1: sp_shot % total shots in plane will change to seq.nshot after disdaqs.
        % ishot
        for nd=1:n_disdaq % n_disdaq will be set to 1 after adc is on.
            % nd
            for isl = 1:seq.nsli
                rf_ex.freqOffset = seq.sl_freqs(isl);
                rf_ex.phaseOffset = rf_ex.phaseOffset+deg2rad(117);
                adc.phaseOffset = adc.phaseOffset+deg2rad(117);
                rf_fatsat.freqOffset = -440;% Hz
                %optional fatsat pulse
                if(adc_flag==0)
                    segmentID = 1;
                else
                    segmentID=2+ishot;
                end

                if(seq.addfatsat)
                    psd.addBlock(rf_fatsat,g_fatsat, mr.makeLabel('SET', 'TRID', segmentID));
                    psd.addBlock(rf_ex,grad_ex,mr.makeDelay(grad_ex.shape_dur+260e-6))
                else
                    %excitation pulse only
                    psd.addBlock(rf_ex,grad_ex,mr.makeDelay(grad_ex.shape_dur+260e-6), mr.makeLabel('SET', 'TRID', segmentID))
                end
                if(adc_flag==0)
                    psd.addBlock(gx_read{1},gy_read{1}) %No ADC
                else
                    psd.addBlock(gx_read{ishot},gy_read{ishot},adc);
                end

                psd.addBlock(g_spoilx,g_spoilz);

                %block for making the TR correct.
                psd.addBlock(mr.makeDelay(delayTR));
            end

        end
        if(adc_flag==0)% After the disdaq loop ends once, set adc to 1.
            adc_flag=1;
            n_disdaq=1;
            sp_shot=seq.nshot;
            break
        end
    end
end
save('seq.mat','seq')

%% Write to .seq file
[ok, error_report]=psd.checkTiming;
psd.plot('timerange',[10*seq.TR/1000 11*seq.TR/1000]);
psd.write(['/home/dinankg/code/psd/pulseq_psd/2d_spiral_gre/',...
    '2D_spiralGRE.seq']);
%     [kfa,ta,kf]=psd.calculateKspacePP();

%% Optional: Save the kspace trajectory:
if(save_ktraj)
dt = adc.dwell;%sec
ktraj = cat(3,cumtrapz(gx),cumtrapz(gy))*dt; %1/m % klen x nshots x 2
ktraj = permute(ktraj,[1,3,2]);% klen x 2 x nshots
n_adc = floor(adc.numSamples/4)*4;
ktraj=ktraj(1:n_adc,:,:); % Cropping only the acquired part
kspace = ktraj.*(2*pi*seq.FOV/seq.npix);
    sp_number=0
    while true
        sp_number = sp_number + 1
        sp_name = ['spiral_ktraj',num2str(sp_number),'.mat'];

        if ~isfile(sp_name)
            sprintf('Saving kspace traj as spiral_ktraj%1d...',sp_number)
            save(sp_name,'kspace','seq')
            break
        end
    end

end
%% Converting to GE readable file
if(strcmp(scanner,'GE'))
    ACTUAL_SCAN = 1
    seq.scanner = 'inside';
    %% Parameters to convert to tar file
    %bins maxview first then goes to next slice index
    sysGE = toppe.systemspecs('maxGrad', 10, ... % G/cm
        'maxSlew', 20, ... % G/cm/ms
        'maxRF', 0.2, ... % Gauss. Must be >= peak RF in sequence
        'rfDeadTime',100,...
        'adcDeadTime', 0, ... % us. Half of 40us since applied both before + after ADC window.
        'psd_rf_wait', 148, ... % RF/gradient delay (us)
        'psd_grd_wait', 156); % ADC/gradient delay (us)

    %% Convert to GE

    seq2ge('2D_spiralGRE.seq',sysGE,'2D_spiralGRE.tar');

    system('tar -xvf 2D_spiralGRE.tar');
    %% Plotting
    figure; toppe.plotseq(sysGE,'timeRange',[20*seq.TR/1000, 21*seq.TR/1000]); % Plot 1 frame
    % for k=100:15:200
    %     toppe.plotseq(sysGE,'timeRange',[seq.TR_eff*k,seq.TR_eff*(k+1)]./1000);
    %     k+[1,5]
    %     drawnow
    % end
    %%
    if ACTUAL_SCAN
        exp_number = 0;
        while true
            exp_number = exp_number + 1;
            exp_name = ['experiments/',datestr(now,'yyyy-mm-dd_'),seq.scanner,'_exp',num2str(exp_number)];
            if ~isfolder(exp_name)
                mkdir(exp_name)
                break
            end
        end
        [status] = copyfile('2D_spiralGRE.tar', ['./',exp_name,'/2D_spiralGRE.tar']);
        [status] = copyfile('2D_spiralGRE.seq', ['./',exp_name,'/2D_spiralGRE.seq']);
        system(append('tar -xvf ','./',exp_name,'/','2D_spiralGRE.tar',...
            ' -C ',exp_name,'/'));
        toppe.utils.coppe('target','inside','cv',789,'use_pw',true);
        save('variables')
        [status] = copyfile('variables.mat', ['./',exp_name,'/variables.mat']);

    end
    delete('*.mod')
end