%% Recon for 3D GRE coded on pulseq
% Dinank Gupta 
% needs MIRT toolbox
%% Loading Data
clear
dname = './data';
flist = dir([dname,'/Scan*.h5']);
file_name = flist(end).name; %Assuming last file is the one we want
fname = [dname,'/',file_name];

if(strcmp(file_name(1),'P'))
    [dat,hdr] = toppe.utils.loadpfile(fname,[],[],[],'acq_order',1);
    dat = flip(dat,1);
else
[dat,hdr] = toppe.utils.loadsafile(fname,'acq_order',1);
end
save_kdata = 1; % Flag to save Raw Data
save_img = 1; % Flag to take 3D NUFFT and save the image
%% Loading the k-space traj
kfolder = './experiments/';

dt = '2024-07-05';
scanner ='inside';
exp_num = '1';
extra = '';
var_name = [kfolder,dt,'_',scanner,'_exp',exp_num,extra,'/variables.mat'];
scan_variables = load(var_name);
seq_params = scan_variables.seq_params;
%% Ktraj
gamT = 4.2576e7;   % Hz/Tesla
gamG = gamT/1e4;    % Hz/Gauss

gx=scan_variables.gx1;gy=scan_variables.gy1;
gy=gy(1:length(gx),:);
dt = scan_variables.adc.dwell;%sec
ktraj = cat(3,cumtrapz(gx),cumtrapz(gy))*dt*gamG*100; %1/m % klen x nshots x 2
ktraj = permute(ktraj,[1,3,2]);% klen x 2 x nshots
n_adc = floor(scan_variables.adc.numSamples/4)*4;
ktraj=ktraj(1:n_adc,:,:); % Cropping only the acquired part
kspace = ktraj.*(2*pi*seq_params.fov(1)/seq_params.N(1));
%% Sorting the kdata and ktrah according to shots
kspace = squeeze(kspace);ks_rot=zeros(size(kspace,1),3,seq_params.N(3),seq_params.nshot_spiral);
n=0;
kdata=[];
for nt = 1: seq_params.ntp
    for ns = 1:seq_params.nshot_spiral
        for kz = 1:seq_params.N(3)
            kz_tmp = seq_params.zp_scale(kz)*repmat(seq_params.kzmax ,[length(kspace),1]);
            kspace(:,3,ns) = kz_tmp;
            ks_rot(:,:,kz,ns) = kspace(:,:,ns);
                n=n+1;
                kdata(:,:,kz,ns,nt) = dat(:,:,n);
        end
    end
end
kspace = ks_rot;

% kdata size is klen x coils x kz x nshots x nc x ntp
% kspace size is klen x 3 x kz x nshots
% Assumes same ktraj across timepoints
if save_kdata
    save([dname,'/',file_name(1:end-3),'_raw'],"kspace","kdata","seq_params")
end
%% 3D recon
if save_img
    % Recon params
    Nx = seq_params.N(1); Ny = seq_params.N(2); Nz = seq_params.N(3);
    nufft_args = {[Ny,Nx,Nz],[6,6,6],[2*Ny,2*Nx,2*Nz],[Ny/2,Nx/2,Nz/2],'table', 2^11, 'minmax:kb'};
    mask = true(Ny,Nx,Nz); % Mask for support
    L = 6;
    kt = reshape(permute(squeeze(kspace),[1,3,4,2]),[],3); % Combining shots. % size is (klen*shots) x 2
    Gn = Gnufft(mask, [{kt}, nufft_args(:)']);
    dcf = pipe_menon_dcf(Gn);
    % Performing 3D NUFFT
    for tp = 1: size(kdata,5) %timepoints loop
        k = kdata(:,:,:,:,tp); %size is klen x coil x kz x rot
        k = permute(k,[1,3,4,2]); %size is klen x kz x rot x coil
        parfor ic = 1:size(k,4) % coils loop
            d2d = k(:,:,:,ic);   % [ndat nleafs] %size is (klen*shots)
            imtmp = Gn'*(d2d(:).*(dcf(:).^(1)));
            imstmp = reshape(imtmp, Nx,Ny,Nz);
            img(:,:,:,ic,tp) = imstmp;
        end
    end
% image size is x-y-z-coil-tps
    % save([dname,'/',file_name(1:end-3),'_img'],"img","seq_params",'-v7.3')
end