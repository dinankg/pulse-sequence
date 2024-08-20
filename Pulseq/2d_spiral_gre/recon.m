%% Recon for 3D GRE coded on pulseq
% Dinank Gupta 
% needs MIRT toolbox
%% Loading Data
clear
dname = 'data/';
flist = dir([dname,'/Scan_20_shots*.h5']);
file_name = flist(end).name; %Assuming last file is the one we want
fname = [dname,'/',file_name];

if(strcmp(file_name(1),'P'))
    [dat,hdr] = toppe.utils.loadpfile(fname,[],[],[],'acq_order',1);
    dat = flip(dat,1);
else
[dat,hdr] = toppe.utils.loadsafile(fname,'acq_order',1);
end
save_kdata = 0; % Flag to save Raw Data
save_img = 1; % Flag to take 3D NUFFT and save the image
%% Loading the k-space traj
kfolder = '../experiments/';
load('spiral_ktraj2.mat')
% dt = '2024-08-19';
% scanner ='inside';
% exp_num = '5';
% extra = '';
% var_name = [kfolder,dt,'_',scanner,'_exp',exp_num,extra,'/variables.mat'];
% scan_variables = load(var_name);
% seq = scan_variables.seq;
%% Ktraj
gamT = 4.2576e7;   % Hz/Tesla
gamG = gamT/1e4;    % Hz/Gauss

%% Sorting the kdata according to shots
% kspace = squeeze(kspace);ks_rot=zeros(size(kspace,1),3,seq.N(3),seq.nshot_spiral);
n=0;
kdata=[];
for nt = 1: seq.ntp
    for ishot = 1:seq.nshot
        for islice = 1:seq.nsli
                n=n+1;
                kdata(:,:,islice,ishot,nt) = dat(:,:,n);
        end
    end
end
% Assumes same ktraj across timepoints
if save_kdata
    save([dname,'/',file_name(1:end-3),'_raw'],"kspace","kdata","seq")
end
%% 2D recon
    % Recon params
    Nx = seq.npix; Ny = seq.npix;
    nufft_args = {[Ny,Nx],[6,6],[2*Ny,2*Nx],[Ny/2,Nx/2],'table', 2^11, 'minmax:kb'};
    mask = true(Ny,Nx); % Mask for support
    L = 6;
    kt = reshape(permute(squeeze(kspace),[1,3,2]),[],2); % Combining shots. % size is (klen*shots) x 2
    Gn = Gnufft(mask, [{kt}, nufft_args(:)']);
    dcf = pipe_menon_dcf(Gn);
    % Performing 2D NUFFT
    for tp = 1: size(kdata,5) %timepoints loop
        k = kdata(:,:,:,:,tp); %size is klen x coil x kz x rot
        for islice = 1:seq.nsli
        parfor ic = 1:size(k,2) % coils loop
            d2d = k(:,ic,islice,:);   % [ndat nleafs] %size is (klen*shots)
            imtmp = Gn'*(d2d(:).*(dcf(:).^(1)));
            imstmp = reshape(imtmp, Nx,Ny);
            img(:,:,islice,ic,tp) = imstmp;
        end
        end
    end
    figure,im(rot90(rms(img(:,:,:,:,1)),2))

% image size is x-y-z-coil-tps
if save_img
    save([dname,'/',file_name(1:end-3),'_img'],"img","seq",'-v7.3')
end