function bval = calc_bval_trap(G,D,d,rt)
%{
Calculate b-value for a trapezoidal gradient 
G: grad amp (G/cm)
D: spacing between grads (ms)
d: length of 1 grad minus the rise time of 1 of the ramps.(ms)
rt: rise_time of the 1 ramp.(ms)
%}
gambar = 4.257*2*pi; %gamma in kHz/Gauss
bval = gambar^2 * G^2 * ((d^2)*(D-(d/3)) + (rt^3)/30 - d*(rt^2)/6);%ms/cm2
bval = bval*(1e-3/1e2)%s/mm2

b_val_rect = gambar^2 * G^2 * d^2 * (D-d/3);
