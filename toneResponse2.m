function toneResponse2
% function toneResponse
%   Outputs tone stimulus using Psychtoolbox's audio driver, as well as a
%   digital syncing signal. Monitors button press response indicating
%   detection of the tone. Digital IO for response and digital sync are
%   handled by a LabJack U3. Saves timestamped stimulus-response data to
%   file.
%
%   DR 06/2021

% default parameters
N = 10; % number of trials
freq = [1000 2000 3000]; % tone frequency (Hz) *can be more than one*
dur = 0.4; % tone duration (s)
amp = 0.15; % tone amplitude (normalized, 0 to 1)
tint = [5 15]; % inter-stimulus time interval (s)
ftap = .05; % fraction of tone length to taper with Tukey window to prevent clicks (0 = no taper, 1 = max taper)
sdir = 'C:\R01NS113366\ORdata\'; % save directory
subj = 'HUPtest'; % subject name

% data directory
cd(sdir);
if exist(subj,'dir')~=7
    mkdir(subj);
end
cd(subj);

% initialize LabJack
NET.addAssembly('LJUDDotNet'); % make the UD .NET assembly visible in MATLAB
ljObj = LabJack.LabJackUD.LJUD;
[~, ljH] = ljObj.OpenLabJackS('LJ_dtU3','LJ_ctUSB', '0', true, 0); % open LabJack
ljObj.ePutS(ljH, 'LJ_ioPIN_CONFIGURATION_RESET', 0, 0, 0); % reset to factory defaults
ljObj.ePutS(ljH, 'LJ_ioPUT_ANALOG_ENABLE_PORT', 0, bin2dec('0000000000001111'), int32(16)); % configure flexible IO as analog in (channel FIO0-3) and digital in/out (channels FIO4-7 + db15 channels) 
ljObj.eDO(ljH, 4, 0); % set FIO4 to output-low (0V)

% initialize Psych-Audio port
InitializePsychSound(1);
hpa = PsychPortAudio('Open'); % open Psych-Audio port (using all default options)
status = PsychPortAudio('GetStatus', hpa);
fs = status.SampleRate;

% create tones of different frequencies
Nfreq = length(freq);
for ii = 1:Nfreq
    tone(ii,:) = tukeywin(fix(dur*fs),ftap)'.*sin(2*pi*freq(ii)*(0:fix(dur*fs)-1)/fs);
end

%Loop here
exitTask = false;
while exitTask == false
    % create buffer with random intervals of silence before each tone
    buffer = ones(1,2*N);
    ifq = randi(Nfreq,1,N); ieven = 1; % randomize tone frequency on each trial, if more than one specified
    isi = tint(1) + diff(tint)*rand(1,N); iodd = 1;
    for ii = 1:2*N
        if mod(ii,2)
            buffer(ii) = PsychPortAudio('CreateBuffer',[],zeros(2,round((isi(iodd)-dur)*fs)));
            iodd = iodd + 1;
        else
            buffer(ii) = PsychPortAudio('CreateBuffer',[],[tone(ifq(ieven),:); tone(ifq(ieven),:)]);
            ieven = ieven + 1;
        end
    end
    buffer(2*N+1) = PsychPortAudio('CreateBuffer',[],zeros(2,round(5*fs))); % add final silence buffer to allow time (5s) for final response
    
    % create schedule
    PsychPortAudio('UseSchedule',hpa,1,2*N+1);
    for ii = 1:2*N+1
        PsychPortAudio('AddtoSchedule',hpa,buffer(ii),1,0,[]); % audio handle, buffer handle, repititions, start sample, end sample = max
    end
    
    % create data matrix
    dat = zeros(N*10,2); % event x [timestamp, code]; code 0 = start tone, 1 = stop tone, 2 = button press, 3 = button release; assumes no more than 10 events on average per trial
    tstart = datestr(now,'yyyy-mm-dd-HH-MM-SS'); % current date-time for filename
    
    % play tone
    PsychPortAudio('Volume', hpa, amp);
    PsychPortAudio('Start', hpa, [], 0, 1); % start immediately, wait for device to really start
    tic; ctrial = 1;
    cbuff = 0; idat = 1;
    cbutton = 1; % pull-up resistor makes button logic-high when not pressed
    while cbuff < 2*N+1
        
        % send out digital sync signal during tone
        status = PsychPortAudio('GetStatus', hpa);
        buff = status.SchedulePosition;
        if buff > cbuff
            if mod(buff,2) && (buff < 2*N+1) % change to odd buffer = start tone
                ljObj.eDO(ljH, 4, 1); % set FIO4 to output-high (3.3V)
                dat(idat,1) = GetSecs; dat(idat,2) = 0; idat = idat + 1;
            elseif ~mod(buff,2) % change to even buffer = stop tone
                ljObj.eDO(ljH, 4, 0); % set FIO4 to output-low (0V)
                dat(idat,1) = GetSecs; dat(idat,2) = 1; idat = idat + 1;
                elapsedTime = toc;
                disp(['trial = ' num2str(ctrial) ', ISIcommand = ' num2str(round(isi(ctrial)*1000)) 'ms, ISIactual = ' num2str(round(elapsedTime*1000)) 'ms']);
                tic
                ctrial = ctrial + 1;
            end
            cbuff = buff;
        end
        
        % response button edge detection
        [~, button] = ljObj.eDI(ljH, 6, 0); % read FIO6
        if button && ~cbutton
            dat(idat,1) = GetSecs; dat(idat,2) = 3; idat = idat + 1;
            disp('button release');
            cbutton = button;
        elseif ~button && cbutton
            dat(idat,1) = GetSecs; dat(idat,2) = 2; idat = idat + 1;
            disp('button press');
            cbutton = button;
        end
    end
    
    % clean up
    PsychPortAudio('Stop', hpa, 1, 1); % wait until playback finishes
    PsychPortAudio('UseSchedule', hpa, 0); % disable schedule
    PsychPortAudio('DeleteBuffer'); % delete buffer
    PsychPortAudio('Close', hpa);
    
    % save data
    dat(idat:end,:) = [];
    fnm = [subj '-' tstart '.mat'];
    freqs = freq(ifq);
    save(fnm,'dat','N','freqs','isi','dur','amp','tint','ftap','subj','-mat');
end
