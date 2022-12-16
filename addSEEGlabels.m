function addSEEGlabels
% function addSEEGlabels
%   Add brain segmentation labels to each sEEG contact and save in the 2nd
%   column of the channel_labels variable. Uses labels found in ITKSnap CSV
%   file.
%
%   DR 12/2022

% parameters
ddir = '/Users/drew/Box/sEEG'; % data directory
subj = 'HUP225_RID700'; % subject

% load locations
cd(ddir);
locCSV = fullfile(subj,[subj '_ITKSNAP'],'electrodenames_coordinates_native_and_T1.csv'); % CSV file path and name
try
    locations = readtable(locCSV);
catch
    error([locCSV ' not found']);
end

% match formatting with what is in channel_labels
locelec = locations.Var1;
for ii = 1:length(locelec)
    if length(locelec{ii}) == 3
        locelec{ii} = [locelec{ii}(1:2), '0', locelec{ii}(3)];
    end
end

% load channel labels
load(fullfile(subj,[subj '_Induction.mat']),'channel_labels'); % the CSV reconstructions are useful only for induction sessions (after EMU stay)

% add location to each channel label
for ii = 1:length(channel_labels)
    elec = channel_labels{ii,1};
    ind = find(strcmp(locelec,elec));
    if ~isempty(ind)
        channel_labels{ii,2} = locations.Var2{ind};
    elseif any(ismember(elec,{'C3','C4','Cz','CZ','F3','F4','Fz','FZ','P3','P4','Pz','PZ'}))
        channel_labels{ii,2} = 'EEG';
    elseif any(ismember(elec,{'LOC','ROC'}))
        channel_labels{ii,2} = 'EOG';
    else
        channel_labels{ii,2} = '';
    end
end

% append to session .mat file
save(fullfile(subj,[subj '_Induction.mat']),'channel_labels','-append'); % overwrites if had run this function previously
