function seegCF(session,tawake,tanesth)
% function seegCF(session,tawake,tanesth)
%   Distribution of eigenvalue criticiality x frequency obtained from
%   stability analysis of sEEG data using multivariate autoregression
%   modeling: individual session analysis.
%
%   session = data structure for recording session from sEEG project.
%   tawake = [start time, stop time] in minutes for awake state.
%   tanesth = [start time, stop time] in minutes for anesthetized state.
%
%   Produces plot in similar style as Solovey et al 2015 J Neurosci Fig. 8
%
%   DR 06/2022

% parameters
lns = 60; % line noise frequencies (Hz; leave empty [] for no filter)
lpc = []; % filter cutoffs (Hz; leave empty [] for no filter)
car = 0; % common-average re-reference? (0 or 1)
twn = 0.5; % time window size (s)
if nargin~=3, error('input session data structure and start/stop times for two brain states'); end

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

% time bins
is = 1:round(twn*fs):size(cdat,1); % window start times (samples)
it = is/fs/60; % window start times (minutes)
Nwin = length(is)-1;

% frequency bins
if isempty(lpc)
    fbn = 1:100;
else
    fbn = lpc(1):lpc(2);
end

% eigenvalue bins
ebn = 0:.01:1.2;

% complex eigenvalue distribution
E1 = []; F1 = []; E2 = []; F2 = [];
for ii = 1:Nwin
    [~,A,cmarg,~,~,~] = arfit(cdat(is(ii):is(ii+1)-1,:),1,1); % first order multivariate AR model
    [~,~,~,~,~,lambda] = armodeFAST(A,cmarg); % eigenvalues
    cE = abs(lambda); % mode stability (could convert this to damping timescale like Solovey)
    cF = abs(angle(lambda))./(2*pi/fs); % mode frequency (Hz)
    if it(ii)>tawake(1) && (it(ii)+twn/60)<tawake(2)
        E1 = [E1; cE]; 
        F1 = [F1; cF];
    elseif it(ii)>tanesth(1) && (it(ii)+twn/60)<tanesth(2)
        E2 = [E2; cE]; 
        F2 = [F2; cF];
    end
end
% i0 = find(E<.01); disp(['removing ' num2str(length(i0)) ' modes (' num2str(round(length(i0)/length(E)*100)) '%) with criticality ~ 0']); E(i0) = []; F(i0) = [];
N1 = histcounts2(E1,F1,ebn,fbn);
pN1 = log10(N1/length(E1));
N2 = histcounts2(E2,F2,ebn,fbn);
pN2 = log10(N2/length(E2));
dN = N1-N2;
ldN = log10(abs(dN));
ldN(dN==0) = 0;
ldN(dN<0) = -ldN(dN<0);
% pN1 = smooth2a(pN1,1,1);
% pN2 = smooth2a(pN2,1,1);
% ldN = smooth2a(ldN,1,1);

% figure
Nclr = 512; hue1 = 0/3; hue2 = 1/3;
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
subplot(1,3,1), hold on;
colormap(gca,hsv2rgb([hue1*ones(Nclr,1), ones(Nclr,1), linspace(0,1,Nclr)']));
imagesc(fbn,ebn,pN1);
axis tight; axis_square(gca);
set(gca,'Box','off','TickDir','out','YDir','normal','FontSize',14);
xlabel('mode frequency (Hz)'); ylabel('mode criticality'); title('awake'); 
hc = colorbar; hc.Label.String = 'log10 (prob)';
subplot(1,3,2), hold on;
colormap(gca,hsv2rgb([hue2*ones(Nclr,1), ones(Nclr,1), linspace(0,1,Nclr)']));
imagesc(fbn,ebn,pN2);
axis tight; axis_square(gca);
set(gca,'Box','off','TickDir','out','YDir','normal','FontSize',14);
xlabel('mode frequency (Hz)'); ylabel('mode criticality'); title({subj,'anesthetized'},'HorizontalAlignment','center');
hc = colorbar; hc.Label.String = 'log10 (prob)'; 
subplot(1,3,3), hold on;
colormap(gca,hsv2rgb([flipud([hue2*ones(Nclr/2,1), ones(Nclr/2,1), linspace(0,1,Nclr/2)']);[hue1*ones(Nclr/2,1), ones(Nclr/2,1), linspace(0,1,Nclr/2)']]));
imagesc(fbn,ebn,ldN);
axis tight; axis_square(gca);
clim = get(gca,'CLim');
set(gca,'Box','off','TickDir','out','YDir','normal','FontSize',14,'CLim',[-max(abs(clim)), max(abs(clim))]);
xlabel('mode frequency (Hz)'); ylabel('mode criticality'); title('awake - anesthetized');
hc = colorbar; hc.Label.String = 'log10 (\Delta)'; 