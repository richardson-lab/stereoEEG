function seegTF(session)
% function seegTF(session)
%   Time-frequency analysis of sEEG data. Computes both average power and
%   average criticality in each time x frequency bin. session = data
%   structure for recording session from sEEG project.
%
%   DR 03/2022

% parameters
lns = 60; % line noise frequencies (Hz; leave empty [] for no filter)
lpc = [5 400]; % filter cutoffs (Hz; leave empty [] for no filter)
car = 1; % common-average re-reference? (0 or 1)
twn = 0.5; % time window size (s)
Nf = 50; % number of frequencies
smt = [2 10]; % smoothing factor [frequency bins, time bins] (leave empty [] for no smoothing)
wavl = 'morl'; % wavelet

% remove bad channels
iDC = find(strncmp(session.channel_labels(:,1),'DC',2))'; % DC channels (not included in bad channel list)
[~,igood] = setdiff(session.channel_labels(:,1),[session.bad_channels session.channel_labels(iDC)]);
clab = session.channel_labels(igood,1);
cdat = session.data(:,igood);
fs = session.sample_rate;
N = length(igood);
subj = [session.subject ' ' session.type];
clear session

% line noise filter
if ~isempty(lns)
    for ii = 1:N
        ln = mtmlinenoise(cdat(:,ii),3,fs,fs,lns);
        cdat(:,ii) = cdat(:,ii) - ln;
    end
end

% lowpass/bandpass filter
if ~isempty(lpc)
    [b,a] = butter(2,lpc/(fs/2));
    cdat = filtfilt(b,a,cdat);
end

% re-reference
if car
    ind = 1;
    while ind<N
        l = clab{ind};
        lead = l(isstrprop(l,'alpha'));
        ilead = find(strncmp(clab,lead,length(lead)));
        disp(['re-referencing ' lead ' x ' num2str(length(ilead))]);
        cdat(:,ilead) = cdat(:,ilead) - mean(cdat(:,ilead),2); % remove common average activity on each lead
        ind = ilead(end) + 1; % assumes labels are sorted
    end
end

% % diagnostics
% pwelchPlot(cdat,fs,[0.5 500],5,1,500,clab); % PSD of cleaned data (takes a long time)
% timePlot(cdat,fs,[],[],1,[],clab); % scrolling time plot of cleaned data

% time bins
is = 1:round(twn*fs):size(cdat,1); % window start times (samples)
Nwin = length(is)-1;
tcenter = is(1:end-1)/fs/60;

% frequency bins
if isempty(lpc)
    fbn = logspace(log10(1),log10(fs/2),Nf+1); % up to Nyquist frequency
elseif numel(lpc)==1
    fbn = logspace(log10(1),log10(lpc),Nf+1); % up to lowpass cutoff
elseif numel(lpc)==2
    fbn = logspace(log10(lpc(1)),log10(lpc(2)),Nf+1); % within bandpass
end
fcenter = (fbn(1:end-1)+fbn(2:end))/2; % center frequencies

% spectrogram
sgram = zeros(Nf,Nwin,N);
scales = centfrq(wavl)./(fcenter'*(1/fs));
for ii = 1:N
    disp(['wavelets ' num2str(ii)]);
    P = cwt(cdat(:,ii),scales,wavl);
    P = abs(P); % power spectra
    for jj = 1:Nwin
        sgram(:,jj,ii) = mean(P(:,is(jj):is(jj+1)-1),2);
    end
end
sgram = trimmean(sgram,10,3); % average across all electrodes

% criticality-ogram
cgram = zeros(Nf,Nwin); 
for ii = 1:Nwin
    [~,A,cmarg,~,~,~] = arfit(cdat(is(ii):is(ii+1)-1,:),1,1); % first order multivariate AR model
    [~,~,~,~,~,lambda] = armodeFAST(A,cmarg); % eigenvalues
    F = abs(angle(lambda))./(2*pi/fs); % mode frequency (Hz)
    [~,~,ix] = histcounts(F,fbn);
    n = zeros(length(fcenter),1);
    for jj = 1:max(ix)
        ind = find(ix==jj);
        if ~isempty(ind)
            n(jj) = mean(abs(lambda(ind))); % mean criticality value within each time x frequency bin
        end
    end
    cgram(:,ii) = n; 
end

% mean spectra plot
smarg = mean(sgram,2);
cmarg = mean(cgram,2);

% conditioning
cgram = cgram-mean(cgram,2); % changes over time relative to mean
sgram = zscore(sgram,0,2);
if ~isempty(smt) % smooth
    cgram = smooth2a(cgram,smt(1),smt(2)); 
    sgram = smooth2a(sgram,smt(1),smt(2));
end 

% mean spectra plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/6 1/4 1/3 1/2],'Color','w');
plot(fcenter,smarg/max(smarg),'b','LineWidth',3); hold on;
plot(fcenter,cmarg,'r','LineWidth',3);
axis tight; set(gca,'Box','off','TickDir','out','FontSize',14,'XScale','log');
xlabel('frequency (Hz)'); ylabel('\color{blue} relative power  \color{red} criticality index');

% TF plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
subplot(2,1,1);
colormap(jet(256));
imagesc(tcenter,1:size(cgram,1),cgram);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
ytick = get(gca,'YTick'); set(gca,'YTickLabel',round(fcenter(ytick)));
xlabel('time (min)'); ylabel('frequency (Hz)'); title(subj);
hc1 = colorbar('Location','eastoutside');
hc1.Label.String = '\Delta criticality';
subplot(2,1,2);
imagesc(tcenter,1:size(sgram,1),sgram);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
ytick = get(gca,'YTick'); set(gca,'YTickLabel',round(fcenter(ytick)));
xlabel('time (min)'); ylabel('frequency (Hz)');
hc2 = colorbar('Location','eastoutside');
hc2.Label.String = '\Delta power (zscore)';
