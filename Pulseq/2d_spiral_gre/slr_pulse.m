function [rf,grad] =  slr_pulse(slthick,tbw,dur,type,ftype,system)
%Temporary code. Will fill in the details.
%%
% slthick : m slice thickness
%%

slthick =slthick*100; % converting to cm
% Input variables:

maxs = 0.9*system.maxSlew/4258/(1e5); %G/cm/ms # 0.9* to make sure to not reach max slew
maxg = system.maxGrad/4258/100; %G/cm
dt = system.rfRasterTime;

%Making RF
rf = dzrf(2+round(dur/dt),tbw,type,ftype ,0.001,0.001).'; % min phase inv pulse
rf = rf(2:end-1);
bw = tbw/dur; %in Hz
rf_len = length(rf);

%Making sl select grad
g_amp = bw/4258/slthick; % Grad amp in G/cm, Using gamma in Hz/G
g_plat = g_amp * ones(rf_len,1);
g_step = maxs*dt*1e3;
g_ramp = linspace(0,g_amp,ceil(g_amp/g_step)).';
grad = cat(1,g_ramp,g_plat,flip(g_ramp));
rf = cat(1,0*g_ramp,rf,0*g_ramp);

% Balancing Gradient
[~,rf_peak_ind]=max(rf);
grad_rephasor_area = sum(grad(rf_peak_ind:end))*dt; %G/cm*s
grad_rephasor = dotrap(grad_rephasor_area*1e4,maxs*10,maxg*10,dt*1e3)/10;

grad = cat(1,grad,-1*grad_rephasor.');
rf = cat(1,rf,0*grad_rephasor.');

end

function grad = dotrap(area,gslew,gamp,gts)
%dotrap - generates a trap waveform
% units must match
% area = mT/m * ms
% gslew = mT/m/ms
% gamp = mT/m
% gts = ms

% determine if trap or tri
maxrampl = ceil(gamp/gslew/gts);
maxtriarea = maxrampl * gamp * gts;
if (area < maxtriarea) 
    rampl = ceil(sqrt(abs(area)/gslew)/gts);
    shape = [(0:rampl) ((rampl-1):-1:0)];
    grad = shape/sum(shape*gts)*area;
else
    flatl = ceil((abs(area)-maxtriarea)/gamp/gts);
    shape = [(0:maxrampl)/maxrampl ones([1 flatl]) ((maxrampl-1):-1:0)/maxrampl];
    grad = shape/sum(shape*gts)*area;
end
end


