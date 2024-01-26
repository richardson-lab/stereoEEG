function pacPreprocess(subject)
% function pacPreprocess
%   Generates a .csv file of sEEG data, with bad channels removed and line
%   noise filter applied. Also generates a .json file containing channel
%   labels, sample rate, and parameters for running permutation test in
%   python
%
%   AG 01/2024

% parameters
fP = logspace(log10(1),log10(20),21); % phase: frequencies (Hz)
fA = logspace(log10(5),log10(200),21); % amplitude: frequencies (Hz)
delta = pi/8; % phase bin size (rad)
tbin = 120; % time bin length (s)
tstp = 30; % time step size (s)
Nperm = 200; % number of permutations for MI significance test (leave empty or set to zero to not compute MIpval which takes a long time)
tperm = [-60 60]; % interval of random permutation time shift (s)

% % initialize parallel computation pool
% parpool(maxNumCompThreads);

% load data

root = '/users/tnl/Library/CloudStorage/Box-Box/sEEG';

matchedDir = fullfile(root, subject);
if ~isfolder(matchedDir)
    error('No matching subdirectory found.');
end

datfile = dir(fullfile(matchedDir, '*_induction.mat'));
if isempty(datfile)
    error('No matching .mat file found in the directory.');
elseif length(datfile) > 1
    error('Multiple matching .mat files found. Expecting only one.');
end

% cd(datfile(ifile).folder); # Restore this after debugging
load(fullfile(matchedDir,datfile(1).name),'bad_channels', ...
    'channel_labels','data','sample_rate');
disp(['Loaded ' datfile(1).name]);

% remove bad channels and white matter channels
ibad = notSEEGchannels(channel_labels,1,bad_channels);
data(:,ibad) = [];
channel_labels(ibad,:) = [];
Nch = size(data,2);

% line noise filter
for ich = 1:Nch
    ln = mtmlinenoise(data(:,ich),3,sample_rate,sample_rate,60:60:max(fA));
    data(:,ich) = data(:,ich) - ln;
end

% re-reference to common average of each lead
ind = 1;
while ind < Nch
    l = channel_labels{ind,1};
    lead = l(isstrprop(l,'alpha'));
    ilead = find(strncmp(channel_labels,lead,length(lead)));
    if length(ilead) > 1
        disp(['re-referencing ' lead ' x ' num2str(length(ilead))]);
        data(:,ilead) = data(:,ilead) - mean(data(:,ilead),2);
    end
    ind = ilead(end) + 1; % assumes labels are sorted
end  

% save
PACparam.fP = fP;
PACparam.fA = fA;
PACparam.delta = delta;
PACparam.tbin = tbin;
PACparam.tstp = tstp;
PACparam.Nperm = Nperm;
PACparam.tperm = tperm;
PACparam.channel_id_labels = channel_labels(:,1);
PACparam.channel_area_labels = channel_labels(:,2);
PACparam.sample_rate = sample_rate;
json_str = jsonencode(PACparam);

saveDir = fullfile('/Users/tnl/matlab/sEEG_proc',subject);
if ~isfolder(saveDir)
    mkdir(saveDir);
end

json_file = [subject '_param.json'];
fid = fopen(fullfile(saveDir,json_file), 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

csv_file = [subject '_data.csv'];
writematrix(data, fullfile(saveDir,csv_file));
disp([subject ' Data conversion complete'])
% delete(gcp);
