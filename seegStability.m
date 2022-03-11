function seegStability(session)
% function seegStability(session)
%   Stability analysis of sEEG data using multivariate autoregression
%   modeling. session = data structure for recording session from sEEG
%   project.
%
%   DR 02/2022

% parameters
lns = 60; % line noise frequencies (Hz; leave empty [] for no filter)
lpc = [5 400]; % filter cutoffs (Hz; leave empty [] for no filter)
car = 1; % common-average re-reference? (0 or 1)
twn = 0.5; % time window size (s)
lbn = 0.05:0.005:1.02; % eigenvalue bins
lth = 0.5; % eigenvalue thresholds (can list more than one)
smt = [2 10]; % smoothing factor [crit index bins, time bins] (leave empty [] for no smoothing)

% remove bad channels
iDC = find(strncmp(session.channel_labels(:,1),'DC',2))'; % DC channels (not included in bad channel list?)
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

% multivariate AR model fit and eigendecomposition
is = 1:round(twn*fs):size(cdat,1); % window start times (samples)
Nwin = length(is)-1;
lambhist = zeros(length(lbn)-1,Nwin);
for ii = 1:Nwin
    [~,A,C,~,~,~] = arfit(cdat(is(ii):is(ii+1)-1,:),1,1); % first order model
    [~,~,~,~,~,lambda] = armodeFAST(A,C);
    lambhist(:,ii) = histcounts(abs(lambda),lbn)';
end
if ~isempty(smt), lambhist = smooth2a(lambhist,smt(1),smt(2)); end

% fraction above threshold 
lambthresh = zeros(length(lth),Nwin);
for ii = 1:length(lth)
    lambthresh(ii,:) = sum(lambhist(lbn(1:end-1)>lth(ii),:))./sum(lambhist);
    if ~isempty(smt), lambthresh(ii,:) = smooth(lambthresh(ii,:),smt(2)); end
end

% plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/6 1/3 2/3],'Color','w');
h1 = subplot(2+length(lth),1,1:2);
colormap(hot(256));
imagesc(is(1:end-1)/fs/60,lbn,lambhist);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
xlabel('min'); ylabel('criticality index'); title(subj);
hc = colorbar('Location','eastoutside');
hc.Label.String = 'number of modes';
for ii = 1:length(lth)
    h(ii) = subplot(2+length(lth),1,2+ii);
    plot(is(1:end-1)/fs/60,lambthresh(ii,:),'k','LineWidth',2);
    axis tight; set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('min'); ylabel(['modes > ' num2str(lth(ii))]);
    l1 = get(h1,'Position'); l2 = get(h(ii),'Position');
    set(h(ii),'Position',[l2(1:2) l1(3) l2(4)]);
end
