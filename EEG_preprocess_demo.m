% The authors of this code are Adrija Chatterjee and Prof. Pragathi P.
% Balasubramani, Translational Neuroscience and Technology Lab(Transit),
% Department of Cognitive Science, Indian Institute of Technology, Kanpur.

%References:
% 1. https://eeglab.org/

% For this code to function efficiently, it will need the plugins- firfilt,
% ICLabel, clean_rawdata-master.

%FILES TO BE ADDED FOR: Line 15, 22, 155
%please uncomment after adding the path, 

% For any query, please contact: Adrija Chatterjee: adrijac23@iitk.ac.in or
%Prof. Pragathi P. Balasubramani: pbalasub@iitk.ac.in 

%% add the packages required to your path. 
% Packages required are uploaded in the same folder 

%addpath('the path');


%% getting chanlocs(channel locations)
%uploading the setfile to be processed.

% setf= pop_loadset('ADD PATH ');

% channel list 
labels = {'fp2','f4','c4','p4','fp1','f3','c3','p3','f8','t4','t6','o2','f7','t3','t5','o1','cz','pz','fz'};
required_labels = {'fp2','f4','c4','p4','fp1','f3','c3','p3','cz','pz','fz'};
labels_to_remove = {'f8','t4','t6','o2','f7','t3','t5','o1'};


%reading the chanlocs txt file 
chansinfo = readmatrix('/Users/adrijachatterjee/Downloads/codes/chanlocs.txt');
temp = readtable('/Users/adrijachatterjee/Downloads/codes/chanlocs.txt');
chanlabels = table2array(temp(:,1));
chansinfo = chansinfo(:,2:end);
chanlocs_new = [];
nchans = length(labels);
%indx=[]
indices_to_remove = []
for i=1:nchans

    idx = find(strcmpi(chanlabels,labels{i}));
    if any(strcmp(labels_to_remove, labels{i}))
        indices_to_remove = [indices_to_remove i];
    end

    chanlocs_new(:,i).X = chansinfo(idx,1);
    chanlocs_new(:,i).Y = chansinfo(idx,2);
    chanlocs_new(:,i).Z = chansinfo(idx,3);

    [sph_theta, sph_phi, chanlocs_new(:,i).sph_radius] = cart2sph(chanlocs_new(:,i).X,chanlocs_new(:,i).Y,chanlocs_new(:,i).Z);
    chanlocs_new(:,i).sph_phi = rad2deg(sph_phi);
    chanlocs_new(:,i).sph_theta = rad2deg(sph_theta);

    [chanlocs_new(:,i).urchanlocs] = i;
    [~,chanlocs_new(:,i).theta,chanlocs_new(:,i).radius] = sph2topo([chanlocs_new(:,i).urchanlocs,chanlocs_new(:,i).sph_phi,chanlocs_new(:,i).sph_theta]);
    chanlocs_new(:,i).labels = char(chanlabels(idx));

    chanlocs_new(:,i).theta = chanlocs_new(:,i).theta + 90 ; %this was done because the positions were rotated by 90 degrees


end
chanlocs_new_removed = chanlocs_new;
setf_removed = setf;
chanlocs_new_removed(:, indices_to_remove) = []; %channel locations of the required channels.
setf_removed(indices_to_remove, :) = []; %data after removing the extra channels.

%% Loading the .mat file into the EEG struct 

            EEG_raw = pop_importdata('dataformat','matlab','nbchan',0,'data',setf_removed,'srate',250,'pnts',0,'xmin',0);
            
%% Creating the EEG struct
           
            EEG_raw.trials = 1;
            EEG_raw.nbchan = size(EEG_raw.data,1);
            EEG_raw.pnts = size(EEG_raw.data,2);
            EEG_raw.srate = 256; %Hz sampling rate
            EEG_raw.xmin = 0;
            EEG_raw.xmax = size(EEG_raw.data,2)/EEG_raw.srate;
            EEG_raw.times = linspace(EEG_raw.xmin,EEG_raw.xmax,EEG_raw.pnts);
            EEG_raw.etc = [];

            %reading channel locations 

            EEG_raw.chanlocs = chanlocs_new_removed;
%% Resampling the EEG_raw data to match that of EGG sampling rate

            EEG_raw = pop_resample(EEG_raw,250);

%% Trimming some data at the end because it is just flat 

            end_rem_dur = 5; %standard
          
            EEG_pop = pop_select(EEG_raw,"time",[EEG_raw.xmin EEG_raw.xmax-end_rem_dur]);
            fprintf('Data removal done');
 %% Filtering the data (FIR filter)

            EEG_filt = pop_eegfiltnew(EEG_pop, 'locutoff',  1, 'hicutoff',  40, 'filtorder', 9000, 'plotfreqz', 0);
            % this gives better results
 %% to visualise in GUI

%setting the normalised data in eeglab 
eeglab
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_raw, CURRENTSET, 'setname', 'new', 'gui', 'off');
eeglab redraw

%% Remove bad channels using clean_channels function

%not sure if this is needed.. need to plot the scroll data and check 

            % EEG_filt_rmchan = clean_channels(EEG_filt);
            % fprintf('Bad channels removed');

%% ICA

EEG_filt = pop_runica(EEG_filt,'runica');
EEG_filt.icachansind = double(1:EEG_filt.nbchan);
EEG_filt = iclabel(EEG_filt);

fprintf('ICA done\n')

%% Remove bad components and reconstruct data

            th_signal = 0.80; % brain component with less than 5% confidence is removed
            classes = EEG_filt.etc.ic_classification.ICLabel.classes; %classes 
            cls_score = EEG_filt.etc.ic_classification.ICLabel.classifications;%classification scores
            bad_comp = []; %to store bad components

            for cmp = 1:size(cls_score, 1)
                if any(cls_score(cmp, 2:6) > th_signal) || any(cls_score(cmp,1)<0.05)
                    bad_comp = [bad_comp, cmp];
                end
              
            end
            EEG_ica = pop_subcomp(EEG_filt, bad_comp, 0);
    
            fprintf('Bad components removed\n');

plots_topo_rm= pop_topoplot(EEG_ica, 0, 1:size(EEG_ica.icaweights, 1), 'Independent Components', 0, 'electrodes', 'on');
%% only if required 

EEG_ica1= clean_rawdata(EEG_ica,5,[0.25 0.75],0.8,4,5,-1);


 %% Interpolate bad channels

            EEG_interpol = pop_interp(EEG_ica1, EEG_raw.chanlocs, 'spherical');
            fprintf('Bad channels interpolated\n');

 %% Re-reference the data
            EEG_reref = pop_reref(EEG_interpol,[]);
            fprintf('Data rereferenced\n');

%% save the file 

%savedf= pop_saveset(EEG_reref,'demo','ADD PATH'); %ADD PATH WHERE YOU WANT TO SAVE YOUR PROCESSED FILE.

pop_saveset(EEG_reref,'savemode','onefile'); 
fprintf('saved dataset\n');     

%% END