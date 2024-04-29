%{
3D SE MR-ARFI scan with a stack of spirals readout.
%}

% Needs toppe v6:

addpath(genpath('/home/dinankg/tools/toppev6')) %get tv6
%

clear

clc
import toppe.*
import toppe.utils.*

%%

seq.FOV = 40;           % Starting Field of view, x-y, cm
seq.nshot = 8;         % For multishot spiral
seq.npix = 256;          % Matrix size, x-y
seq.und_samp_fact = 1 % Undersampling factor. Should be divisble by nshot for now.

seq.cv_num = str2double(sprintf('%d', 70, seq.FOV)); %enrty file number.

%% Setup sequence parameters using getparams function

% KEY parameters set here

seq.zFOV = 4.5;          % z FOV, cm
seq.nsli = 1;            % Number of slices
seq.zsteps = 15; %
seq.ntp = 3;
seq.spiraldir = 1;     % 1 is spiral out, 2 is spiral in
sp_sl_sc = .7;  % vds spiral slew scaling 0.8 (and presumably higher) could
% crash bc of system hardware limits, 0.7 no warnings or crash

seq.dt = 4e-6;          % TOPPE sampling rate in seconds

seq.fa=90;
seq.slthick = seq.zFOV;     % Slice thickness in cm %for 3D excitation
seq.tbw = 6;           % Time bandwidth for SLR pulse
seq.pulsedur = 4;      % RF pulse duration in ms
seq.rf_sl_sc = 0.76;    % rf slew scale for calling makeslr function


seq.TR=500;%TR
% TE is set to the min TE based on the gradients used.

%% setup toppe system using defaults

% system = toppe.systemspecs();
system = toppe.systemspecs('maxSlew', 20);
ncycles = 0; % balanced sl select.
    system.maxGrad = 10;
    % system.maxGrad = 5;
system.maxRF = 0.25;
%% make tipdown.mod

if seq.fa < 20, type = 'st'; else, type = 'ex'; end
freq_ss=0;
[rf, gex, freq_ss, fnamestem] = toppe.utils.rf.makeslr(seq.fa, seq.slthick, seq.tbw, ...
    seq.pulsedur, ncycles, system,...
    'maxSlewScale',seq.rf_sl_sc,'type',type,'ofname', 'tipdown.mod');
% Append 0s to account for nChop (psd_rf_wait) 25 samples=100us
rf = cat(1,zeros(25,1),rf);
gex = cat(1,zeros(25,1),gex);
seq.gamp90 = max(gex);

toppe.writemod(system, 'rf', rf, 'gz', gex,'nChop',[25,0], 'ofname', 'tipdown.mod');
%% make rephasor.mod
ncycles_180 = 20;
seq.tbw180 = 4;
seq.pulsedur180 = 6;
freq_ss=0;
seq.slfact_180 =.9;
[rf_180, gex_180, freq_180, fnamestem_180] = toppe.utils.rf.makeslr(180,...
    seq.slfact_180*seq.slthick, seq.tbw180, ...
    seq.pulsedur180, ncycles_180, system,...
    'maxSlewScale',seq.rf_sl_sc*0.8,'type','se','ofname', 'tmp.mod',...
    'discardPrephaser',1,'ftype','ls');
delete('tmp.mod')
if(max(rf_180)>system.maxRF)
    error('RF peak too large');
end

seq.gamp180 = max(gex_180(length(gex_180)/2-50:length(gex_180)/2+50));
toppe.writemod(system, 'rf', rf_180, 'gz', gex_180, 'ofname', 'rephasor.mod');


%% Adding bipolar ARFI grads to this:
dwamp = 4; %G/cm
dwlen = 3 * 1e-3; %sec
dw_dir=[0 1 0]; % direction to play yhe ARFI gradient.
dw_gradient = toppe.utils.trapwave2(dwamp*dwlen, dwamp, 15, seq.dt*1e3);
dw_gx=dw_dir(1)*dw_gradient';
dw_gy=dw_dir(2)*dw_gradient';
dw_gz=dw_dir(3)*dw_gradient';

% Delay between the 2 gradients to allow ultrasound to rise to max
seq.dw_grad_delay = 1.5; %ms
grad_delay_samp = seq.dw_grad_delay /seq.dt/1e3;
seq.num_dw_grad = 2; % For OGSE type grads. 2 sets it to bipolar on both sides of 180
% Flip and append
dw_gxall = dw_gx; dw_gyall = dw_gy; dw_gzall = dw_gz;
for g =1:seq.num_dw_grad - 1
    dw_gxall = [dw_gxall;zeros(grad_delay_samp,1);dw_gx(2:end)*(-1)^g];
    dw_gyall = [dw_gyall;zeros(grad_delay_samp,1);dw_gy(2:end)*(-1)^g];
    dw_gzall = [dw_gzall;zeros(grad_delay_samp,1);dw_gz(2:end)*(-1)^g];
end
% Time to send the trigger (us):
quick_trig=0;
if(seq.num_dw_grad==1 && quick_trig==1)
    t_trigout = 200;
elseif(seq.num_dw_grad==1)
    t_trigout = 1000;
else
    t_trigout = length(dw_gradient)* seq.dt*1e6+100;
end
seq.grad_dir = [1 1 -1]; 
% ^If you want to flip the polarity of the ARFI grad at any timepoint,set that to -1 here.
toppe.writemod(system,'gx',dw_gxall,'gy',dw_gyall,'gz',dw_gzall,'ofname','arfigrad.mod');
%% Making list of triggers
seq.triglist = ones(1,seq.ntp);
% Turning off for first time point
seq.triglist(1,[1]) = 0;

%% Generate balanced VD spiral

seq.oversamp = 600;

[gx,gy,t,npts] = makevdspiral2(seq.FOV,seq.npix,seq.nshot,seq.oversamp,5,sp_sl_sc*system.maxSlew*10);
paramsint16(2) = npts;
gz=0*gx;
gx=gx(:,1);gy=gy(:,1);gz=gz(:,1);
if(seq.spiraldir==2)
    fprintf('Not tested with spiral in')
    gx=flip(gx);gy=flip(gy);gz=flip(gz);
end
% List of rotations for spirals:
shotrot = linspace(0,2*pi-(2*pi)/seq.nshot*seq.und_samp_fact,seq.nshot/seq.und_samp_fact);
seq.rot_list =cat(2,repmat([0 0 1],seq.nshot/seq.und_samp_fact,1),shotrot.');
seq.rotmat = axang2rotm(seq.rot_list);

%% Make steps in z
gambar = 4257;                           % gamma/2pi in Hz/Gauss
dkz = 1/seq.zFOV; %1/cm
res = seq.zFOV/seq.zsteps;%cm
seq.kzmax = 1/res;%1/cm
Gz_area = seq.kzmax/gambar; %G sec/cm grad area to reach kz max

% Grads will go from -Gz_area/2 to Gz_area/2
zp = (-(seq.zsteps-1)/2:(seq.zsteps-1)/2)/((seq.zsteps-1)/2)/2; %zpe scaling for each TR
%List of gradient area:
gz_trap_area_list = -1*Gz_area; %G sec/cm
%First trap which will be the longest one.
z_enc_trap = toppe.utils.trapwave2(gz_trap_area_list(1),6, system.maxSlew-5,system.raster*1e-3).';
seq.zp_scale = zp; %Scaling of the kz traps
% Appending z enc to spirals.
gz=cat(1,z_enc_trap,gz,-1*z_enc_trap);
gx=cat(1,0*z_enc_trap,gx,0*z_enc_trap);
gy=cat(1,0*z_enc_trap,gy,0*z_enc_trap);
paramsint16(14) = length(z_enc_trap);
%NOTE than total number of samples to skip will be paramsint16(14)+paramsint16(11)
% And the spiral in samples are paramsint16(2) for recon.
toppe.writemod(system,'gx',gx,'gy',gy,'gz',gz,'ofname','readout.mod','hdrints',paramsint16);
%% Make a table of rotation matrices and kz-encode scaling.
GR = 1.618;
seq.doGA=0; % Flag to turn GA on/off
GA_tp_n = randi(1e2,seq.ntp,1); % GA differences across time.
seq.rot_angle_shot=[];seq.z_scale_shot=[];% shots x tp
for nt=1:seq.ntp
    seq.n_shots_per_vol=0;
    for nspiral = 1:(seq.nshot/seq.und_samp_fact)
        n=0; %Resetting the golden angle rotations for the same kz platters,
        %just adding rotations
        for nkz = 1:seq.zsteps
            seq.n_shots_per_vol=seq.n_shots_per_vol+1;
            n=n+1;
            seq.rot_angle_shot(seq.n_shots_per_vol,nt) = mod(seq.doGA*(360-(360/GR))*(n+GA_tp_n(nt))+...
                rad2deg(shotrot(nspiral)),360);
            seq.z_scale_shot(seq.n_shots_per_vol,nt) = seq.zp_scale(nkz);
            % seq.rot_angle_shot(seq.n_shots_per_vol,nt)
        end
    end
end
seq.rot_angle_shot = deg2rad(seq.rot_angle_shot);
seq.rot_list_shot =cat(2,repmat([0 0 1],seq.n_shots_per_vol*seq.ntp,1),seq.rot_angle_shot(:));
seq.rotmat_shot = axang2rotm(seq.rot_list_shot);
seq.rot_list_shot = reshape(seq.rot_list_shot,[],seq.ntp,4);%shots x tp x 4
seq.rotmat_shot = reshape(seq.rotmat_shot,3,3,[],seq.ntp);% 3 x3 x shots x tp
%% Crusher at the end of spiral.
gspoil = toppe.utils.makecrusher(4, seq.slthick/seq.zsteps, system, 0, 4,system.maxGrad);
toppe.writemod(system,'gx',gspoil,'gy',gspoil,'gz',gspoil,'ofname','crusher.mod','hdrints',paramsint16);
%% Delays to get the right timings:
fprintf('Still recommend seeing the waveforms on the oscilloscope to ensure that the delays make a spin-echo.\n')
% Left of 180 should match right of 180 (and probably remove extra time
% that DAQ core takes to begin).
n_rf_left = length(rf) - (seq.pulsedur/2)/seq.dt/1e3;
seq.delay_90_180 = length(z_enc_trap)*seq.dt*1e3; % ms Can add delays needed for fmaps too
seq.delay_180_ro = seq.dt*n_rf_left*1e3 - length(z_enc_trap)*seq.dt*1e3;
time_current = (length([rf;dw_gx;dw_gxall;rf_180;gx(:,1);gspoil])+n_rf_left)*seq.dt*1e3;
time_current =time_current+seq.delay_90_180 +seq.delay_180_ro;
time_add = seq.TR - time_current;
seq.delay_TR = time_add;

%% Making modules and cores.txt
% Entries are tab-separated.
% Calculating the min durations for each mod:
tipdown_dur = length(rf)*seq.dt*1e6;
rephasor_dur = length(rf_180)*seq.dt*1e6;
readout_dur = length(gx)*seq.dt*1e6;
crusher_dur = length(gspoil)*seq.dt*1e6;
arfi_dur = length(dw_gxall)*seq.dt*1e6;
modFileText = ['' ...
    'Total number of unique cores\n' ...
    '5\n' ...
    'fname	duration(us)	hasRF?	hasDAQ? hastrigout?\n' ...
    'tipdown.mod	',num2str(tipdown_dur),'	1	0   -1\n'...
    'rephasor.mod	',num2str(rephasor_dur),'	1	0   -1\n'...
    'readout.mod	',num2str(readout_dur),'	0	1   -1\n' ...
    'crusher.mod	',num2str(crusher_dur),'	0	0   -1\n'...
    'arfigrad.mod ',num2str(arfi_dur),'   0   0   ' num2str(t_trigout)];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);


coreFileText{1} = [1,5,0,2,5];%tip-grad-<delay>-180-grad
coreFileText{2} = [0,3];%<delay>-kz-spiral-kz
coreFileText{3} = [4,0];%crush-<delay>

writecoresfile(coreFileText);


%% Setup the sequence timing


writeloop_3d_arfi(freq_ss, seq, system);


%% Scan time:
printf(['Time per image volume: ',num2str(seq.TR*seq.nshot*seq.zsteps/1e3),' sec'])
%% Play sequence out virtually
figure;

for j=1:11
    toppe.plotseq(1+9*(j-1),9*j,system,'doDisplay','true','gmax',dwamp);
    drawnow
end
%%

% figure,
% % m=toppe.utils.rf.slicesim([0 0 1],rf_b,gz_b,4e-3,-3:0.05:3,1000,100);

