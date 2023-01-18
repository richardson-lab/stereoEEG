function ibad = badChannels(channel_labels,wm,varargin)
% function ibad = badChannels(channel_labels,wm,bad_channels)
%   Returns indices of 'bad' sEEG channels. It will look through
%   channel_labels and find all the channels that typically are not used in
%   our sEEG analyses (e.g. scalp EEG channels) and return indices for
%   these. If second input is set to 1, will also exclude channels with
%   white matter labels. If third input is given, it will also exclude
%   channels listed in bad_channels cell array.
%   
%   channel_labels (ch x 2 cell array)
%   wm = 0 or 1 (include or exclude white matter channels?)
%   bad_channels = cell array of previously identified bad channel labels
%
%   DR 01/2023

ibad = [];
for ii = 1:length(channel_labels)
    elec = channel_labels{ii,1};
    labl = channel_labels{ii,2};
    if any(ismember(elec,{'C3','C4','Cz','CZ','Fp1','Fp2','F3','F4','F7','F8','Fz','FZ','P3','P4','Pz','PZ','O1','O2','T3','T4','T5','T6'})) % exclude EEG
        ibad = [ibad ii];
    elseif any(ismember(elec,{'LOC','ROC'})) % exclude EOG
        ibad = [ibad ii];
    elseif strncmp(elec,'EMG',3) % exclude EMG
        ibad = [ibad ii];
    elseif strncmp(elec,'EKG',3) % exclude EKG
        ibad = [ibad ii];
    elseif strncmp(elec,'SMELL',4) % exclude SMELL signal
        ibad = [ibad ii];
    elseif strncmp(elec,'RESP',4) % exclude RESP signal
        ibad = [ibad ii];
    elseif strncmp(elec,'DC',2) % exclude DC (analog input channels)
        ibad = [ibad ii];
    elseif isempty(labl) % exclude channels without label (usually outside of brain)
        ibad = [ibad ii];
    elseif wm && contains(labl,'White Matter')
        ibad = [ibad ii];
    end
end
if nargin==3
    bad_channels = varargin{1};
    bad = zeros(length(channel_labels),1);
    for ii = 1:length(bad_channels)
        bad = bad + strcmp(channel_labels(:,1),bad_channels{ii});
    end
    ibad = [ibad, find(bad)'];
end
ibad = sort(unique(ibad));