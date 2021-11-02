% Downloads data from iEEG.org and saves in a structure
% Use iEEG portal to select a start and end time for saved segment
% Optional: Enter noisy channels noticed on iEEG into bad_channels array
% For description of saved variables, see last segment of code

%% Variables
% Machine/user variables
ieeg_pw_file = '/Users/drew/Box/Matlab/ieeg/agr_ieeglogin.bin'; % directory of password file created by IEEGSession.createPwdFile
ieeg_user = '********';
data_dir = '/Users/drew/Data/R01NS113366'; % save location for downloaded data

% Session variables
subject = 'HUP223'; % subject HUP###
session_type = 'Emergence'; % Induction or Emergence
session_title = 'HUP223_OR_implant_emergence';
start_seconds = 447; % manually selected start time
end_seconds = 3050; % manually selected end time
bad_channels = {'LA07','LD12'}; % Manually enter bad channels

%% Data download
% Start iEEG Session
ieeg_session = IEEGSession(session_title, ieeg_user, ieeg_pw_file);

% Build raw_data matrix (samples x chanels)
[~,num_channels] = size(ieeg_session.data.rawChannels);
start_useconds = start_seconds*1000000;
end_useconds = end_seconds*1000000;
raw_data = getvalues(ieeg_session.data,start_useconds,end_useconds-start_useconds,1); % getvalues(data, start microsec, duration, channel)
for i=2:num_channels
    raw_data(:,i) = getvalues(ieeg_session.data,start_useconds,end_useconds-start_useconds,i);
end
t = [0:size(raw_data,1)-1] / ieeg_session.data.sampleRate;

% Read Annotations, offset to match extracted data time
ieeg_annotations = getEvents(ieeg_session.data.annLayer, 0);
annotations = {ieeg_annotations(1,1).description,(ieeg_annotations(1,1).start- start_useconds)/1000000};
for i=2:length(ieeg_annotations)
    annotations = [annotations;{ieeg_annotations(1,i).description,(ieeg_annotations(1,i).start - start_useconds)/1000000}];
end

% Read targets from spreadsheet
if isfile(fullfile(data_dir,subject,strcat(subject, '_', session_type, '_channels.xlsx')))
channel_labels = ieeg_session.data.channelLabels;
[~,~,sheet] = xlsread(fullfile(data_dir,subject,strcat(subject,  '_', session_type, '_channels.xlsx')));
else
    disp('No channel spreadsheet found. lead_targets will be left blank')
    sheet = {};
end

%% Concatenate toneResponse data
if strcmp(session_type, 'Induction')
    tone_files = dir(fullfile(data_dir,subject,'toneResponse_Induction','*.mat'));
    load(fullfile(tone_files(1).folder,tone_files(1).name));
    tones.dat = dat;
    tones.N = N;
    tones.freqs = freqs;
    tones.isi = isi;
    tones.dur = dur;
    tones.amp = amp;
    tones.tint = tint;
    tones.ftap = ftap;
    tones.subj = subj;
    for i=2:length(tone_files)
        load(fullfile(tone_files(i).folder,tone_files(i).name));
        tones.dat = [tones.dat; dat];
        tones.N = tones.N + N;
        tones.freqs = [tones.freqs, freqs];
        tones.isi = [tones.isi, isi];
    end
else
    tones = 'no toneResponse data collected';
end

%% Build and save session structure
session.subject = subject;                                  % HUP### assigned to subject
session.type = session_type;                                % Emergence or Induction
session.sample_rate = ieeg_session.data.sampleRate;         % sample rate
session.start_us = start_useconds;                          % offset in milliseconds from start of recording
session.end_us = end_useconds;                              % milliseconds at end of downloaded segment
session.channel_labels = ieeg_session.data.channelLabels;   % labels for each lead/channel
session.data = raw_data;                                    % raw eeg data in microvolts
session.t = t;                                              % time in seconds for extracted data
session.annotations = annotations;                          % annotations typed during session
session.lead_targets = sheet;                               % anatomical target of each lead
session.bad_channels = bad_channels;                        % manually identified noise channels
session.tones = tones;                                      % concatenated toneResponse data

save(fullfile(data_dir,subject,strcat(subject,'_', session_type, '.mat')),'session','-V7.3');
