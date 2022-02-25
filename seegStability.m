function seegStability(session)
% function seegStability(session)
%   Stability analysis of sEEG data using multivariate autoregression
%   modeling. session = data structure for recording session from sEEG
%   project.
%
%   DR 02/2022

% parameters
lpc = 100; % low-pass filter cutoff (Hz; leave empty [] for no filter)
lns = 60; % line noise frequencies (Hz; leave empty [] for no filter)
car = 1; % common-average re-reference? (0 or 1)
twn = 0.5; % time window size (s)
lbn = 0.8:0.005:1.02; % eigenvalue bins
lth = 0.98; % eigenvalue threshold
smt = 2; % smoothing factor (leave empty [] for no smoothing)

% remove bad channels
iDC = find(strncmp(session.channel_labels(:,1),'DC',2))'; % DC channels (not included in bad channel list?)
[~,igood] = setdiff(session.channel_labels(:,1),[session.bad_channels session.channel_labels(iDC)]);
clab = session.channel_labels(igood,1);
cdat = session.data(:,igood);
fs = session.sample_rate;
N = length(igood);
subj = [session.subject ' ' session.type];
clear session

% low-pass filter
if ~isempty(lpc)
    [b,a] = butter(2,lpc/(fs/2));
    cdat = filtfilt(b,a,cdat);
end

% line noise filter
if ~isempty(lns)
    for ii = 1:N
        ln = mtmlinenoise(cdat(:,ii),3,fs,fs,lns);
        cdat(:,ii) = cdat(:,ii) - ln;
    end
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
lambhist = zeros(length(lbn)-1,N);
for ii = 1:Nwin
    [~,A,C,~,~,~] = arfit(cdat(is(ii):is(ii+1)-1,:),1,1); % first order model
    [~,~,~,~,~,lambda] = armodeFAST(A,C);
    lambhist(:,ii) = histcounts(abs(lambda),lbn)';
end
if ~isempty(smt), lambhist = smooth2a(lambhist,smt,1); end

% fraction above threshold 
lambthresh = sum(lambhist(lbn(1:end-1)>lth,:))./sum(lambhist);
if ~isempty(smt), lambthresh = smooth(lambthresh,2*smt); end

% plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
h1 = subplot(3,1,1:2);
colormap(hot(256));
imagesc(is(1:end-1)/fs/60,lbn,lambhist);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
xlabel('min'); ylabel('criticality index'); title(subj);
hc = colorbar('Location','eastoutside');
hc.Label.String = 'number of modes';
h2 = subplot(3,1,3);
plot(is(1:end-1)/fs/60,lambthresh,'k','LineWidth',2);
axis tight; set(gca,'Box','off','TickDir','out','FontSize',14);
xlabel('min'); ylabel(['modes > ' num2str(lth)]);
l1 = get(h1,'Position'); l2 = get(h2,'Position');
set(h2,'Position',[l2(1:2) l1(3) l2(4)]);
