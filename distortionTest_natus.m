function distortionTest_natus
% function distortionTest_natus
%   Analyze frequency response data to assess any signal amplitude/phase
%   distortions due to filtering by the data acquisition system. Assumes
%   inputing constant-amplitude, variable-frequency sine wave into an
%   amplifier channel (called 'response' below) and a sync signal with the
%   zero crossings of that sine wave into a digital input channel. 
%
%   DR 03/2022

% parameters
ddir = '/Users/drew/Box/sEEG/Frequency Response Test/'; % data dir
dfile = '20220325_4096.edf'; % data file (either .edf or .mat)
chsync = 'DC1'; % name of sync channel
chresp = 'RI2';%'CH001';%'RI2'; % name of response data channel
fs = 4096; % sampling rate (Hz)
ttl = 0.1; % threshold for TTL signal (V)
ampi = 5000; % input sine wave amplitude (uVpp)
np = []; % number of poles in transfer function estimate
freqs = [0.05, 0.1:0.1:0.9, 2:200]; % frequencies to analyze (Hz)

% load data
cd(ddir);
if contains(dfile,'edf') 
    [hdr, record] = edfread(dfile);
    try
        isync = find(strcmp(hdr.label,chsync));
        sync = record(isync,:)/1e6; % V
    catch
        error(['no channel named ' chsync]);
    end
    try
        iresp = find(strcmp(hdr.label,chresp));
        resp = -1*record(iresp,:); % uV **invert since otherwise phase lags didn't make sense (i.e. 180deg out of phase in passband)**
    catch
        error(['no channel named ' chresp]);
    end
    clear record
elseif contains(dfile,'mat')
    load(dfile,'session');
    try
        isync = find(strcmp(session.channel_labels,chsync));
        sync = session.data(:,isync)/1e6; % V
    catch
        error(['no channel named ' chsync]);
    end
    try
        iresp = find(strcmp(session.channel_labels,chresp));
        resp = -1*session.data(:,iresp); % uV **invert since otherwise phase lags didn't make sense (i.e. 180deg out of phase in passband)**
    catch
        error(['no channel named ' chresp]);
    end
end
if contains(dfile,'20220322_1024')
    istop = (21*60+19)*fs; % stop sample--after this, recording was bad
    resp(istop:end) = [];
    sync(istop:end) = [];
end
inan = find(isnan(resp));
resp(inan,:) = [];
sync(inan,:) = [];

% frequency and gain from hilbert transform of response
hresp = hilbert(resp);
amp = abs(hresp);
phase = angle(hresp);
uphase = unwrap(phase);
ipeaks = round(interp1(uphase,1:length(phase),0:2*pi:floor(uphase(end)/(2*pi))*2*pi))'; % peak of each cycle (sample)
ipeaks = ipeaks(2:end-1); 
freq = fs./diff(ipeaks); % frequency of each cycle (Hz)
ind = find(freq<0.1); freq(ind) = round(freq(ind),2);
ind = find(freq>=0.1 & freq<1); freq(ind) = round(freq(ind),1);
ind = find(freq>=1); freq(ind) = round(freq(ind));
g = 2*amp(ipeaks(1:end-1))'/ampi; % gain of each cycle

% phase response inferred from sync signal
icyc = find(sync(1:end-1)<ttl & sync(2:end)>=ttl); % beginning of each sine cycle (sample; could be delayed by up to 1/fs sec due to sampling rate)
d = NaN*ones(size(g));
for ii = 1:length(ipeaks)-1
    [~,jj] = min(abs(icyc-ipeaks(ii))); % closest sync zero crossing to response peak
    ipki = icyc(jj)+fs/freq(ii)/4; % peak of input signal (occuring at cycle period/4) in samples
    if abs(ipki-ipeaks(ii))>(fs/freq(ii))/4 % closest input peak more than 90 deg away from response peak
        continue;
    end
    d(ii) = (ipki-ipeaks(ii))*360/(fs/freq(ii)); % deg
end

% statistics across cycles at same frequency
gain = NaN*ones(length(freqs),2);
delay = NaN*ones(length(freqs),2);
for ii = 1:length(freqs)
    jj = find(freq==freqs(ii)); 
    if isempty(jj)
        continue;
    end
    glim = prctile(g(jj),[5 95]);
    kk = find(g(jj)>glim(1) & g(jj)<glim(2));
    gain(ii,:) = [mean(g(jj(kk))) std(g(jj(kk)))]; % trim mean and trim std
    mu = angle(nansum(exp(1i*d(jj)*pi/180)))*180/pi; % delay is a circular variable
    delay(ii,:) = [mu, 0]; % std for circular variable?
end
inan = find(isnan(delay(:,1)));
gain(inan,:) = [];
delay(inan,:) = [];
freqs(inan) = [];

% plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
subplot(2,1,1), errorbar(freqs,gain(:,1),gain(:,2),'ko','LineStyle','none','MarkerFaceColor','k');
axis tight; set(gca,'Box','off','YGrid','on','TickDir','out','FontSize',14,'XScale','log');
ylabel('gain'); xlabel('frequency (Hz)'); title(dfile,'Interpreter','none');
subplot(2,1,2), errorbar(freqs,delay(:,1),delay(:,2),'ko','LineStyle','none','MarkerFaceColor','k');
axis tight; set(gca,'Box','off','YGrid','on','TickDir','out','FontSize',14,'XScale','log');
ylabel('phase (deg)'); xlabel('frequency (Hz)');

% estimate transfer function
if ~isempty(np)
    response = gain(:,1).*exp(1i*delay(:,1)*pi/180);
    gfr = idfrd(response,freqs*2*pi,1/fs);
    sys = tfest(gfr,np)
    [mag,ph] = bode(sys,freqs*2*pi);
    subplot(2,1,1), hold on; plot(freqs,squeeze(mag),'b');
    subplot(2,1,2), hold on; plot(freqs,squeeze(ph),'b');
end
