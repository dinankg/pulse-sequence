 function M = bloch2(Mi, Beff, gamma, dt, nstep)
%function M = bloch(Mi, Beff, gamma, dt, nstep)
%	Evolve magnetization field using time-discretized Bloch equation
%	Input
%		Mi	[dim,3]		initial X,Y,Z magnetization
%					Note: normalized to magnitude <= Mo=1
%		Beff	[dim,3]		effective X,Y,Z applied magnetic field
%					(Tesla), for rotating frame (no Bo)
%		gamma	[dim]		gyromagnetic ratio (MHz/Tesla)
%		dt	scalar		time interval (msec)
%		nstep	scalar		number of time steps to do (default=1)
%	Output
%		M	[dim,3]		final X,Y,Z magnetization
%
%	The "dim" dimensions can be anything provided all are consistent.
%
%	Magnetization is expressed in a frame of reference rotating at the
%	Larmor frequency of water having gyromagnetic ratio of 42.57 MHz/Tesla,
%	assuming a static z-directed magnetic field of 1.5 Tesla.
%
%	Rewritten 10/19/97 for Matlab5 by J. Fessler
%	Adapted 11/17/94 by M. Lubinski
%	Original 11/5/86 bloch.f by cj hardy (ge-crd)

% Check function call
if nargin < 4, error('Not enough input arguments'), end
if nargin < 5, nstep = 1; end

ndim = ndims(Mi);

% Here are some debugging checks.  Turn on to check, off for speed.
if any(size(Mi) ~= size(Beff)), error('M Beff size mismatch'), end
dims = size(Mi);
dimg = size(gamma);
if any(dims(1:end-1) ~= dimg(1:end-1)), error('M gamma size mismatch'), end

%
% reshape to be [n,3] - would be nice to make this more elegant.  ideas?
% maybe use subsref?  the problem is want Beff(:,:,...,:,3) - how?
%
%allbutlast = cells(1, ndim-1);
%for ii=1:ndim-1
%	allbutlast{ii} = 1:dims(ii);
%end
%

Mi = reshape(Mi, prod(dims)/3, 3);
Beff = reshape(Beff, prod(dims)/3, 3);
gamma = gamma(:);

%
%	Adjust effective Bz to account for differences in gamma from H2O
%	question: why aren't all three components scaled?
%
gamma_H2O = 42.57;		% MHz/Tesla
B0 = 3;%1.5;			% Tesla
Beff(:,3) = Beff(:,3) + (gamma/gamma_H2O-1) * B0;

%
%	Adjust "effective B" (really omega) for units to save flops
%
B = dt * 1e3 * 2*pi * Beff;

%
%	Compute sines & cosines of field angles:
%	Theta = angle w.r.t positive z axis
%	Phi   = angle w.r.t positive x axis
%	Psi   = angle w.r.t transformed positive x axis
%
Bmag = sqrt(sum(B.^2,2));		% Magnitude of applied field
Btrans = sqrt(B(:,1).^2 + B(:,2).^2);	% Magnitude of transverse applied field

ct = ones(size(B,1),1);
good = Bmag ~= 0;
if any(good)
	ct(good) = B(good,3) ./ Bmag(good);	% cos(theta)
end
st = sqrt(1 - ct.^2);				% sin(theta) > 0

cphi = ones(size(B,1),1);
good = Btrans ~= 0;
if any(good)
	cphi(good) = B(good,1) ./ Btrans(good);	% cos(phi)
end
sphi = sqrt(1 - cphi.^2) .* sign(B(:,2));	% sin(phi)

cpsi = cos(gamma .* Bmag);			% cos(psi)
spsi = sin(gamma .* Bmag);			% sin(psi)


%
%	Evolve
%
M = Mi;
for ii = 1:nstep

	if 0
		% this was not accurate enough unfortunately, for moderate dt
		M = M + cross(M, B, 1);
	else
		Mx0 = M(:,1);
		My0 = M(:,2);
		Mz0 = M(:,3);

M(:,1) = cphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
+ spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
- sphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
+ cpsi.*(cphi.*My0-sphi.*Mx0));
M(:,2) = sphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
+ spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
+ cphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
+ cpsi.*(cphi.*My0-sphi.*Mx0));
M(:,3) = ct.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0)) ...
- st.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
+ spsi.*(cphi.*My0-sphi.*Mx0));

	end

end

	M = reshape(M, dims);	% return output shape to be same as input
