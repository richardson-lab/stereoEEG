function syncToneResponse
% function syncToneResponse
%   Convert the the tone and response times (see toneResponse.m) from
%   behavioral laptop timestamps to sEEG Natus timestamps and save as a new
%   field in the 'session' structure.
%
%   DR 11/2021

% parameters
ddir = '/Users/drew/Box/sEEG/HUP224'; % data directory
sess = 'HUP224_Induction.mat'; % session data

% check for tone-response data
cd(ddir);
load(sess,'-mat','session');
if ~isfield(session,'tones')
    error('no toneResponse behavioral data detected in session structure');
end
iDC01 = find(strcmp('DC01',session.channel_labels(:,1)),1); % assumes sync signal recorded on DC01 of Natus
if isempty(iDC01)
    error('no DC01 channel found containing sync signal');
end
if isfield(session.tones,'ndat')
    warning('already ran this function');
end

% synchronize
iOnB = find(session.tones.dat(:,2)==0); % tone onset: behavioral laptop
iOnN = find(session.data(1:end-1,iDC01)<1e6 & session.data(2:end,iDC01)>=1e6); % tone onset: natus (assumes 1 V threshold captures binary sync signal state; data in uV)
tOnB = session.tones.dat(iOnB,1); % time of tone onset (s)
tOnN = session.t(iOnN)';
dtOnB = diff(tOnB); % intervals between consecutive tone onsets (s)
dtOnN = diff(tOnN);
[r,lags] = xcorr(dtOnB,dtOnN); % cross-correlation between intervals
[~,imx] = max(r);
lmx = lags(imx);
isyncB = lmx+1:length(dtOnB); % synchronized tone-onset indices
isyncN = 1:min(length(dtOnN),length(isyncB));
[b,~,~,~,stats] = regress(tOnN(isyncN),[ones(length(isyncB),1), tOnB(isyncB)]); % linear mapping from laptop timestamps to Natus timestamps

% verify synchronization
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/8 1/3 3/4 1/3],'Color','w');
subplot(1,3,1); hold on;
plot(lags,r,'k');
axis tight; axis square; set(gca,'Box','off','TickDir','out','XGrid','on');
text(lmx,r(imx),{[],num2str(lmx)},'HorizontalAlignment','center','VerticalAlignment','bottom');
xlabel('lag (intervals)'); ylabel('cross correlation');
subplot(1,3,2); hold on;
plot(dtOnB(isyncB),dtOnN(isyncN),'ko');
ierr = find(abs(dtOnB(isyncB)-dtOnN(isyncN))>0.01); % interval error > 10 ms
if ~isempty(ierr)
    plot(dtOnB(isyncB(ierr)),dtOnN(isyncN(ierr)),'ro','MarkerFaceColor','r');
    title('syncronization not successful');
else
    title('syncronization successful');
end
line([min(dtOnB),max(dtOnB)],[min(dtOnB),max(dtOnB)],'Color','k');
axis square; set(gca,'Box','off','TickDir','out');
xlabel('tone intervals: laptop (s)'); ylabel('tone intervals: natus (s)');
subplot(1,3,3); hold on;
plot(tOnB(isyncB),tOnN(isyncN),'ko');
line(get(gca,'XLim'),b(1)+b(2)*get(gca,'XLim'),'Color','r');
title(['t_{natus} = ' num2str(b(1)) ' + ' num2str(b(2)) ' * t_{laptop} (r^2 = ' num2str(stats(1)) ')']);
axis square; set(gca,'Box','off','TickDir','out');
xlabel('tone onset: laptop (s)'); ylabel('tone onset: natus (s)');
if ~isempty(ierr), error('syncronization not successful'); end

% put all behavioral events on Natus timestamp and save
ind = find(session.tones.dat(:,1)>=tOnB(isyncB(1)) & session.tones.dat(:,1)<tOnB(isyncB(end))); % maybe don't need this, but to be safe consider only events within time range of sync
session.tones.ndat = session.tones.dat(ind,:);
session.tones.ndat(:,1) = b(1)+b(2)*session.tones.ndat(:,1); % add as a new field to the session.tones sub-structure
code = {'tone on';'tone off';'button press';'button release'}; % also add as a new expanded annotations field to session structure
session.annotations2 = [session.annotations; [code(session.tones.ndat(:,2)+1), num2cell(session.tones.ndat(:,1))]];
save(sess,'session','-V7.3');
