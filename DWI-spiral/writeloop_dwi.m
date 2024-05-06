function writeloop_dwi(freq_ss,seq,system)

% Inputs:
%   freq_ss     slice select frequency scale
%   seq         sequence information



%% Calculate frequency offsets for each slice

% calculate bandwidth of slice select pulse in Hz
sli_bw = seq.tbw/(1e-3*seq.pulsedur);
sli_bw180 = seq.tbw180/(1e-3*seq.pulsedur180);
sli_bw180_shift = system.gamma*seq.gamp180*seq.slthick*(1-seq.slfact_180)/2;
%% Adding spacing between adjacent slices:
slispacing_bw90 = round(system.gamma*seq.gamp90*seq.slspacing);
slispacing_bw180 = round(system.gamma*seq.gamp180*seq.slspacing);

%^ find out what needs to change here to get the correct slice plan
% find offset of 1st slice
freq_ss_start = (freq_ss - ( (sli_bw+slispacing_bw90) * (seq.nsli-1)));
freq_ss180_start = (freq_ss - ((sli_bw180+slispacing_bw180) * (seq.nsli-1)));

% create vector with all slice offsets, rounding to nearest Hz
freq_offset_all = -freq_ss_start/2 + round(freq_ss_start:(sli_bw+slispacing_bw90):freq_ss);
freq_offset_all180 = round(-freq_ss180_start/2 + round(freq_ss180_start:(sli_bw180+slispacing_bw180):freq_ss));

if numel(freq_offset_all) ~= seq.nsli
    error('Error calculating slice offsets, check freq_ss parameter.')
end


%% Loop through all timepoints, slices, and Nc indices



ver = 5;
toppe.write2loop('setup',system,'version',ver);
k=0;
seq.disdaqtp=1;
itp=0;

for bv = 1:length(seq.bval_bvec_list.bval)
    grad_amp = seq.bval_bvec_list.G_ratio_list(bv)*seq.bval_bvec_list.bvec(:,bv);

    itp=itp+1;
    freq_offset = freq_offset_all;
    freq_offset180 = freq_offset_all180;
    
    for ishot = seq.disdaqtp:seq.nshot % start from 0 for disdaq
        %         ishot
        for n = 1:seq.nsli
            k=k+1;
            
            if(ishot<1)
                is=1;
                dab='off';
            else
                is=ishot;
                dab='on';
            end
            %once ishot==1, change the disdaqtp to 1 so that next tp
            %onwards it begins as 1.
            if(ishot==1 && seq.disdaqtp==0)
                seq.disdaqtp=1;
            end
            
            
            if(seq.dofatsat)
                % volume fat sat RF pulse & crushers
                toppe.write2loop('fatsat.mod',system,...
                    'RFoffset',-440,'echo',is,'slice',...
                    n,'view',itp,'version',ver);
            end
            %90
            toppe.write2loop('tipdown.mod',system,'RFoffset',...
                freq_offset(n),'echo',is,'slice',n,...
                'view',itp,'version',ver);
            
            % dw grad
            toppe.write2loop('dwgrad.mod',system,...
                'echo',is,'slice',n,'view',itp,...
                'version',ver,'trigout',0,...
                'Gamplitude',grad_amp);
            %180
            toppe.write2loop('rephasor.mod',system,'RFoffset',...
                freq_offset180(n),'echo',is,'slice',n,...
                'view',itp,'version',ver);
            
            % dw grad
            toppe.write2loop('dwgrad.mod',system,...
                'echo',is,'slice',n,'view',itp,...
                'version',ver,'textra',seq.t_extra_arfi2,...
                'Gamplitude',grad_amp);
            
            % spiral readout
            toppe.write2loop('readout.mod',system,...
                'echo',is,'slice',n,'view',itp,...
                'waveform',is,'version',ver,'dabmode',dab);
            
            toppe.write2loop('crusher.mod',system,...
                'echo',is,'slice',n,'view',itp,...
                'version',ver,'textra', seq.Tpad,...
                'trigout', 0);
            
            
        end
    end
end


toppe.write2loop('finish',system);


% Create 'sequence stamp' file for TOPPE.
% This file is listed in line 6 of toppe0.entry
toppe.preflightcheck(['toppe',num2str(seq.cv_num),'.entry'], 'seqstamp.txt', system);

end


