function seegStability(session)
% function seegStability(session)
%   Stability analysis of sEEG data using multivariate autoregression
%   modeling: individual session analysis. 
%   
%   Produces plots in same style as Solovey et al 2015 J Neurosci
%   (Figs. 2, 3, 4). session = data structure for recording session from
%   sEEG project.
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
sdir = '/Users/drew/Data/'; % save directory (leave empty [] to not save anything)

% remove bad channels
iDC = find(strncmp(session.channel_labels(:,1),'DC',2))'; % DC channels (not included in bad channel list?)
[~,igood] = setdiff(session.channel_labels(:,1),[session.bad_channels session.channel_labels(iDC)]);
clab = session.channel_labels(igood,1);
cdat = session.data(:,igood);
fs = session.sample_rate;
N = length(igood);
subj = [session.subject '_' session.type];
annotations = session.annotations;
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
tcenter = is(1:end-1)/fs+twn/2; % s

% multivariate AR model fit and eigendecomposition
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
imagesc(tcenter/60,lbn,lambhist);
set(gca,'Box','off','TickDir','out','FontSize',14,'YDir','normal');
xlabel('min'); ylabel('criticality index'); title(subj,'Interpreter','none');
hc = colorbar('Location','eastoutside');
hc.Label.String = 'number of modes'; h = zeros(1,length(lth));
for ii = 1:length(lth)
    h(ii) = subplot(2+length(lth),1,2+ii);
    plot(tcenter/60,lambthresh(ii,:),'k','LineWidth',2);
    axis tight; set(gca,'Box','off','TickDir','out','FontSize',14);
    xlabel('min'); ylabel(['modes > ' num2str(lth(ii))]);
    l1 = get(h1,'Position'); l2 = get(h(ii),'Position');
    set(h(ii),'Position',[l2(1:2) l1(3) l2(4)]);
end

% save data for across-subject analyses (only saves first threshold data; overwrites any prior saves for this subject)
if ~isempty(sdir)
    try
       load([sdir 'seegStability.mat'],'stability');
       ind = find(strcmp({stability.subj},subj), 1);
       if isempty(ind)
           ind = length(stability)+1;
       end
    catch
        ind = 1;
    end
    stability(ind).subj = subj;
    stability(ind).lns = lns;
    stability(ind).lpc = lpc;
    stability(ind).car = car;
    stability(ind).twn = twn;
    stability(ind).lth = lth(1);
    stability(ind).smt = smt;
    stability(ind).anno = annotations;
    stability(ind).t = tcenter;
    stability(ind).lambthresh = lambthresh(1,:);
    try
        save([sdir 'seegStability.mat'],'stability','-append');
    catch
        save([sdir 'seegStability.mat'],'stability');
    end
end
