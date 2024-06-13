%%
%%
fprintf('This is for TOPPE v5.')
clear
clc
import toppe.*
import toppe.utils.*


%%
seq.FOV = 40;           % Starting Field of view, x-y, cm
seq.nshot = 8;         % For multishot spiral
seq.npix = 256;          % Matrix size, x-y

seq.cv_num = str2double(sprintf('%d', 76, seq.FOV)); %enrty file number.

%% Setup sequence parameters using getparams function

% KEY parameters set here

seq.zFOV = 4.5;          % z FOV, cm
seq.nsli = 15;            % Number of slices

seq.slspacing = .1; % Spacing between adjacent slices, cm
% Note that if slspacing>0 then zFOV is larger than prescribed.
seq.ntp = 1;
seq.spiraldir = 1;     % 1 is spiral out, 2 is spiral in
sp_sl_sc = .7;  % vds spiral slew scaling 0.8 (and presumably higher) will
% crash bc of system hardware limits, 0.7 no warnings or crash


%% Key timing variables
seq.dt = 4e-6;          % TOPPE sampling rate in seconds
seq.TR = 1500;            % TR in ms
seq.TR_eff = seq.TR/seq.nsli; %For interleaving slices

%% RF settings
seq.slthick = seq.zFOV/seq.nsli;     % Slice thickness in cm
seq.dofatsat = 0;
%% Throw away this many scans before steady state
seq.disdaqs = 0;        % discarded timepoints
%% System parameters
seq.maxs_xy = 18;       % max slew in x-y dims, in Gauss/cm/ms
%SlicePlanning inputs
seq.sliceOffset0 = (seq.zFOV/2) - (seq.slthick/2); % assume slices are centered at magnet center
% Scan plane prescription (rotation and slice offset)
% scan at iso-center, non-oblique
seq.rotmat = eye(3);
seq.topSliceOffset = seq.sliceOffset0;

seq.fa=90;
% seq.slthick = seq.zFOV;     % Slice thickness in cm %for 3D excitation
seq.topSliceOffset = seq.topSliceOffset; % in cm


%% Create modules.txt, with entries for water RF, and readout

% Entries are tab-separated.
modFileText = ['' ...
    'Total number of unique cores\n' ...
    '7\n' ...
    'fname	duration(us)	hasRF?	hasDAQ? hastrigout?\n' ...
    'readout.mod	0	0	1   -1\n' ...
    'tipdown.mod	0	1	0   -1\n'...
    'fatsat.mod	0	1	0   -1\n'...
    'rephasor.mod	0	1	0   -1\n'...
    'crusher.mod	0	0	0   -1\n'...
    'dwgrad.mod 0   0   0   -1\n'];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);


%% setup toppe system using defaults

% system = toppe.systemspecs();
system = toppe.systemspecs('maxGrad',5,'maxSlew', 14, 'slewUnit', 'Gauss/cm/ms','maxRF',0.24);
ncycles = 0; % balanced sl select.
% To get timing right:
system.timessi = 200;
system.myrfdel = 148;
system.daqdel = 156;

%% make fatsat.mod 
% fat excitation/saturation settings

seq.fa_fat = 90;            % flip angle in degrees
seq.slthick_fat = seq.zFOV;%slthick; % in cm
seq.bw_desired_fat = 440;   % fat frequency offset in Hz
seq.tbw_fat = 2;            % time bandwidth product
seq.ncycles_fatsat = 20;

% set pulse dur in ms
pulsedur_fat = 1000 * (seq.tbw_fat/seq.bw_desired_fat);

% create RF excitation at fat frequency
[rf_fat, gex2, freq2, fnamestem2] = toppe.utils.rf.makeslr(seq.fa_fat, ...
    seq.slthick_fat, seq.tbw_fat, pulsedur_fat, eps, system, 'type', 'ex',...
     'ofname','tmp.mod'); 
delete('tmp.mod');

% create crushers to go with fat sat
gx_crush0 = toppe.utils.makecrusher(seq.ncycles_fatsat,seq.slthick_fat, system, 0, system.maxSlew-8, system.maxGrad);
% toppe.utils.makecrusher(8, seq.slthick, system, 0, system.maxSlew-8, system.maxGrad); 
gx_crush = [zeros(numel(rf_fat),1); gx_crush0];
gy_crush = gx_crush;
gz_crush = gx_crush;

% pad RF to have time for crushers
rf_fat = [rf_fat; zeros(numel(gx_crush0),1)];
% Adding 250 0s to fat sat to incorporate the trigger.
% rf_fat = [zeros(1000,1);rf_fat];
% gx_crush = [zeros(1000,1);gx_crush];

toppe.writemod(system,'rf', rf_fat,'gx',gx_crush,'ofname','fatsat.mod')

%% make tipdown.mod
seq.pulsedur = 6;      % RF pulse duration in ms
seq.tbw = 6;           % Time bandwidth for SLR pulse
seq.rf_sl_sc = 0.76;    % rf slew scale for calling makeslr function

if seq.fa < 20, type = 'st'; else, type = 'ex'; end
freq_ss=0;
[rf, gex, freq_ss, fnamestem] = toppe.utils.rf.makeslr(seq.fa, seq.slthick, seq.tbw, ...
    seq.pulsedur, ncycles, system, 'sliceOffset', seq.topSliceOffset,...
    'maxSlewScale',seq.rf_sl_sc,'type',type,'ofname', 'tmp.mod','discardPrephaser',1);
delete('tmp.mod')

seq.gamp90 = max(gex);

toppe.writemod(system, 'rf', rf, 'gz', gex, 'ofname', 'tipdown.mod');

%% make rephasor.mod
ncycles_180 = 4;
seq.tbw180 = 6;
seq.pulsedur180 = 6;
freq_ss=0;
seq.slfact_180 =1;
[rf_180, gex_180, freq_180, fnamestem_180] = toppe.utils.rf.makeslr(180,...
            seq.slfact_180*seq.slthick, seq.tbw180, ...
            seq.pulsedur180, ncycles_180, system, 'sliceOffset', seq.topSliceOffset,...
            'maxSlewScale',seq.rf_sl_sc,'type','se','ofname', 'tmp.mod',...
            'discardPrephaser',1,'ftype','ls');
delete('tmp.mod')

seq.gamp180 = max(gex_180(length(gex_180)/2-50:length(gex_180)/2+50));
toppe.writemod(system, 'rf', rf_180, 'gz', gex_180, 'ofname', 'rephasor.mod');

%% Adding bipolar ARFI grads to this:
dwamp = 5; %G/cm
dwlen = 17 * 1e-3; %sec
dw_gradient = toppe.utils.trapwave2(dwamp*dwlen, dwamp, 8, seq.dt*1e3);
dw_1=dw_gradient';

% Delay between the 2 gradients to allow ultrasound to rise to max
seq.dw_grad_delay = 0; %ms
grad_delay_samp = seq.dw_grad_delay /seq.dt/1e3;
seq.num_dw_grad = 1; % For OGSE type DWI.
% Flip and append
dw_all = dw_1;
for g =1:seq.num_dw_grad - 1
dw_all = [dw_all;zeros(grad_delay_samp,1);dw_1(2:end)*(-1)^g];
end
% Time to send the trigger (us):
t_trigout = length(dw_gradient)* seq.dt*1e6;

toppe.writemod(system,'gx',dw_all,'gy',dw_all,'gz',dw_all,'ofname','dwgrad.mod');

% Adding extra time delay before and after the 180
%% Making list of triggers
seq.triglist = ones(1,seq.ntp);
% Turning off for first and last 2 tps
seq.triglist(1,1:2) = 0;
seq.triglist(1,end-1:end) = 0;

%% New readmod with the correct delay for trig

% Entries are tab-separated.

modFileText = ['' ...
    'Total number of unique cores\n' ...
    '6\n' ...
    'fname	duration(us)	hasRF?	hasDAQ? hastrigout?\n' ...
    'readout.mod	0	0	1   -1\n'...
    'fatsat.mod	0	1	0   -1\n' ...
    'tipdown.mod	0	1	0   -1\n'...
    'rephasor.mod	0	1	0   -1\n'...
    'crusher.mod	0	0	0   -1\n'...
    'dwgrad.mod 0   0   0   200'];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

%% Time to append to second arfigrad to make the pulseseq Spin Echo:

% Left of 180 should match right of 180 (and probably remove extra time
% that DAQ core takes to begin).
n_rf_left = length(rf) - (seq.pulsedur/2)/seq.dt/1e3 - 8; % Time from midpoint of rf to its end.
% n_rf_left = (seq.pulsedur/2)/seq.dt/1e3 - 8;
seq.t_extra_arfi2 = seq.dt*n_rf_left*1e3;
% Since DAQ is delayed due to all instructions, subtracting that time from
% here:
seq.t_extra_arfi2 =seq.t_extra_arfi2 - 0.220;

%% Calculate the b value and set the diffusion directions
ramp_time = find(diff(dw_all)==0,1)*seq.dt*1e3;
d = length(dw_all)*seq.dt*1e3 - ramp_time;
D = (length(dw_all)+ length(rf_180))*seq.dt*1e3;
bval = calc_bval_trap(dwamp,D,d,ramp_time);
%List of b values
seq.bval_list = [0,500,bval];
%Scaling the gradient amp to get the b-values in the list.
seq.G_ratio_list = sqrt((dwamp^2)/bval*seq.bval_list)/dwamp;

% setting diffusion directions -right now just x,y,z
I=eye(3);
seq.bval_bvec_list.bval = [0,seq.bval_list(2)*ones(1,3),seq.bval_list(3)*ones(1,3)];
seq.bval_bvec_list.G_ratio_list = [0,seq.G_ratio_list(2)*ones(1,3),seq.G_ratio_list(3)*ones(1,3)];
seq.bval_bvec_list.bvec = [I(:,1),repmat(I,1,2)];

assert(length(seq.bval_bvec_list.bval)==length(seq.bval_bvec_list.bvec),'bval and bvec should be same size')

%% Generate balanced VD spiral
seq.oversamp = 600;

[gx_grad,gy_grad,t,npts] = makevdspiral2(seq.FOV,seq.npix,seq.nshot,seq.oversamp,5,sp_sl_sc*system.maxSlew*10);
paramsint16(2) = npts;% p 
toppe.writemod(system,'gx',gx_grad,'gy',gy_grad,'gz',0*gy_grad,'ofname','readout.mod','hdrints',paramsint16);

%% Crusher at the end of spiral.
gspoil = toppe.utils.makecrusher(8, seq.slthick, system, 0, ...
                    system.maxSlew-10, system.maxGrad-4); 

toppe.writemod(system,'gx',gspoil,'gy',gspoil,'gz',gspoil,...
    'ofname','crusher.mod');

%% Time to append to end of 1 TR to make interleave TR correct.
time_current = (length([rf;dw_1;dw_1;rf_180;gx_grad(:,1);gspoil])+n_rf_left)*seq.dt*1e3;
time_add = seq.TR_eff - time_current;
if(time_add<0)
    error('Increase TR  or decrease slice number to allow interleaving')
end
%% Call writeloop_fmn to create full sequencesystem.slewUnitsystem.slewUnit

seq.Tpad = time_add;

writeloop_dwi(freq_ss, seq, system);

%% Play sequence out virtually
figure(1);pause(1)
toppe.playseq(6+seq.dofatsat,system,'drawpause',1,'tpause',.1,...
    'nTRskip',15,'gmax',system.maxGrad)
% 
figure('Renderer', 'painters', 'Position', [10 10 1200 800])
toppe.plotseq(1,6+seq.dofatsat,system,'doDisplay','true','gmax',10);

arrayfun(@(x) grid(x,'on'), findobj(gcf,'Type','axes'))
arrayfun(@(x) grid(x,'minor'), findobj(gcf,'Type','axes'))

set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',10, ...
'FontWeight','Bold', 'LineWidth', 1,'layer','top');%grid on
