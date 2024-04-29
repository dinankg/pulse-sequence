function writeloop_3d_arfi(freq_ss,seq,system)
% Make the scanloop file for 3D ARFI scan
% Dinank Gupta, April 2024
% University of Michigan

% Inputs:
%   freq_ss     slice select frequency scale
%   seq         sequence information
%   system      system information

% Outputs:
% scanloop.txt, seqstamp.txt and toppe#.entry


%% Calculate frequency offsets for each slice

% calculate bandwidth of slice select pulse in Hz
sli_bw = seq.tbw/(1e-3*seq.pulsedur);
sli_bw180 = seq.tbw180/(1e-3*seq.pulsedur180);
%% Adding spacing between adjacent slices:
% This is not really used for 1 slice/slab scans.
slispacing_bw90 = round(system.gamma*seq.gamp90);
slispacing_bw180 = round(system.gamma*seq.gamp180);

% find offset of 1st slice
freq_ss_start = (freq_ss - ( (sli_bw+slispacing_bw90) * (seq.nsli-1)));
freq_ss180_start = (freq_ss - ((sli_bw180+slispacing_bw180) * (seq.nsli-1)));

% create vector with all slice offsets, rounding to nearest Hz
freq_offset_all = -freq_ss_start/2 + round(freq_ss_start:(sli_bw+slispacing_bw90):freq_ss);
freq_offset_all180 = round(-freq_ss180_start/2 + round(freq_ss180_start:(sli_bw180+slispacing_bw180):freq_ss));

if numel(freq_offset_all) ~= seq.nsli
    error('Error calculating slice offsets, check freq_ss parameter.')
end

%% Loop through all timepoints

ver = 6;
dab='off';
seq.n_disdaq=10;
adc_flag=0;
toppe.write2loop('setup',system,'version',ver);
for itp = 0:seq.ntp
    n=0;

    for ishot = 1:seq.n_shots_per_vol
        for nd = 1:seq.n_disdaq % n_disdaq will be set to 1 during adc

            kz = seq.z_scale_shot(ishot,max(1,itp));
            rotmatx = seq.rotmat_shot(:,:,ishot,max(1,itp));
            n=n+1;

            toppe.write2loop('tipdown.mod',system,'RFoffset',freq_offset_all, ...
                'echo',1,'slice',n,'view',1,...
                'version',ver, 'core',1);

            % dw grad
            toppe.write2loop('arfigrad.mod',system,...
                'echo',1,'slice',n,'view',1,...
                'version',ver,'trigout', seq.triglist(max(1,itp)),...
                'Gamplitude',[1;1;1]*seq.grad_dir(max(1,itp)), 'core',1);
            % Delay module to get the correct SE timing
            toppe.write2loop('delay',system,'textra',seq.delay_90_180,'core',1)

            %180 mod
            toppe.write2loop('rephasor.mod',system,'RFoffset',freq_offset_all180, ...
                'echo',ishot,'slice',n,'view',1,...
                'version',ver,'core',1);

            toppe.write2loop('arfigrad.mod',system,...
                'echo',1,'slice',n,'view',1,...
                'version',ver,'trigout', 0,...
                'Gamplitude',[1;1;1]*seq.grad_dir(max(1,itp)), 'core',1);

            % Delay module to get the correct SE timing
            toppe.write2loop('delay',system,'textra',seq.delay_180_ro,'core',2)
            % spiral readout
            toppe.write2loop('readout.mod',system,...
                'echo',1,'slice',n,'view',max(1,itp),...
                'rotmat',rotmatx,'version',ver,...
                'Gamplitude',[1;1;adc_flag*kz],...
                'dabmode',dab, 'core',2);

            toppe.write2loop('crusher.mod',system,...
                'echo',1,'slice',n,'view',1,...
                'version',ver, 'core',3);

            % Delay module to get the correct TR
            toppe.write2loop('delay',system,'textra',seq.delay_TR,'core',3)

        end
        if(adc_flag==0)% After the disdaq loop ends once, set adc to 1.
            adc_flag=1;
            seq.n_disdaq=1;
            dab='on';
            break
        end

    end

end

toppe.write2loop('finish',system);

toppe.preflightcheck(['toppe',num2str(seq.cv_num),'.entry'], 'seqstamp.txt', system);
copyfile(['toppe',num2str(seq.cv_num),'.entry'],'toppeN.entry')
end