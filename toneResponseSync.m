function toneResponseSync
% function toneResponseSync
%   Convert the the tone and response times (see toneResponse.m) from
%   behavioral laptop timestamps to sEEG Natus timestamps and append to the
%   annotations variable.
%
%   DR 11/2021

% parameters
ddir = '/Users/drew/Box/sEEG/HUP230_RID808/HUP230_Induction'; % data directory
sess = 'HUP230_Induction.mat'; % session data
synch = 'DC01'; % name of channel receiving sync pulses in toneResponse.m

% check for tone-response data
cd(ddir);
load(sess,'-mat');
if ~isstruct(tones)
    error('no toneResponse behavioral data detected in session data');
end
isync = find(strcmp(synch,channel_labels(:,1)),1);
if isempty(isync)
    error(['no ' synch ' channel found containing sync signal']);
end
if isfield(tones,'ndat')
    warning('already ran this function');
end

% synchronize
iOnB = find(tones.dat(:,2)==0); % tone onset: behavioral laptop
iOnN = find(data(1:end-1,isync)<1e6 & data(2:end,isync)>=1e6); % tone onset: natus (assumes 1 V threshold captures binary sync signal state; data in uV)
tOnB = tones.dat(iOnB,1); % time of tone onset (s)
tOnN = t(iOnN)';
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
ierr = find(abs(dtOnB(isyncB)-dtOnN(isyncN))>0.1); % interval error > 100 ms
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
ind = find(tones.dat(:,1)>=tOnB(isyncB(1)) & tones.dat(:,1)<tOnB(isyncB(end))); % maybe don't need this, but to be safe consider only events within time range of sync
tones.ndat = tones.dat(ind,:);
tones.ndat(:,1) = b(1)+b(2)*tones.ndat(:,1); % add as a new field to the session.tones sub-structure
code = {'tone on';'tone off';'button press';'button release'}; % also add as a new expanded annotations
for ii = 1:length(code) % remove prior expanded annotations if running this function again
    annotations(strcmp(annotations(:,1),code{ii}),:) = [];
end
annotations = [annotations; [code(tones.ndat(:,2)+1), num2cell(tones.ndat(:,1))]];
annotations = sortrows(annotations,2);
save(sess,'tones','annotations','-append');
