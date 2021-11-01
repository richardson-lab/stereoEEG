function timePlot(dat,fs,varargin)
% function timePlot(dat,fs,fA,fB,dns,events,chn)
%   Scrolling time-domain plot.
%
%   DR 05/2018

% parameters
fA = [0.3 300]; % filter A: lowpass/bandpass filter cutoff(s) or 'mua'
fB = []; % filter B: [], lowpass/bandpass filter cutoff(s), or 'mua'
dns = 1; % downsample factor (<2 = don't downsample)
events = []; % events specified either as a structure (events.samp (Nx[sample,eventcode]), events.name) or cell array {event, time of occurance in seconds}
chn = []; % channel names
if nargin>=3, fA = varargin{1}; end
if nargin>=4, fB = varargin{2}; end
if nargin>=5, dns = varargin{3}; end
if nargin>=6, events = varargin{4}; end
if nargin>=7, chn = varargin{5}; end

% extract A and B features
DATa = dat;
if ~isempty(fB), DATb = DATa; else DATb = []; end
clear dat;
for ii = 1:2
    if ii==1, cf = fA; cDAT = DATa; else cf = fB; cDAT = DATb; end
    if ~isempty(cf)
        if isnumeric(cf) % lowpass/bandpass filter
            [b,a] = butter(2,cf/(fs/2));
            cDAT = filtfilt(b,a,cDAT);
        elseif strcmp(cf,'mua') % multiunit activity
            [b,a] = butter(2,[300 6000]/(fs/2)); % spike band
            cDAT = filtfilt(b,a,cDAT);
            for ich = 1:size(cDAT,2) % clip signal at 2*std
                cdat = cDAT(:,ich);
                clip = 2*median(abs(cdat)/0.6745);
                cdat(cdat>clip) = clip;
                cdat(cdat<-clip) = -clip;
                cDAT(:,ich) = cdat;
            end
            [b,a] = butter(2,50/(fs/2)); % low pass the rectified signal
            cDAT = filtfilt(b,a,abs(cDAT));
        end
    end
    
    % downsample
    if dns>1
        cDAT = downsample(cDAT,dns); % downsample
    end
    if ii==1, DATa = cDAT; else DATb = cDAT; end
end
if dns>1, fs = round(fs/dns); end
clear cDAT

% interactive plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4],'Color','w');%,'Renderer','Painters'); % ADD THIS TO BE ABLE EXPORT FIGURE AS VECTOR GRAPHICS
ha = axes('Box','off','Units','normalized','TickDir','out','FontSize',12); xlabel('sec'); ylabel('channel'); hold on;
set(ha,'YLim',[0 size(DATa,2)+1],'YTick',1:size(DATa,2));
if ~isempty(chn) && length(chn)==size(DATa,2), set(ha,'YTickLabel',chn); end
pos = get(ha,'Position');
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025]);
hc = uicontrol('Style','listbox','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.1, 0.05, 0.1],'String',{'60','20','10','5','2','1','0.5'},'Value',2,'Max',1,'Min',1,'Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4), 0.05, 0.02],'String','sec','BackgroundColor',get(gcf,'Color'));
he(1) = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.15, 0.05, 0.025],'String','1000','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3)+0.05, pos(2)+pos(4)-0.15, 0.02, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'));
he(2) = uicontrol('Style','edit','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.2, 0.05, 0.025],'String','100','ForegroundColor','b','Callback',@PlotCallback);
uicontrol('Style','text','Units','normalized','Position',[pos(1)+pos(3)+0.05, pos(2)+pos(4)-0.2, 0.02, 0.02],'String','uV','BackgroundColor',get(gcf,'Color'),'ForegroundColor','b');
hb(1) = uicontrol('Style','checkbox','Units','normalized','Position',[pos(1)+pos(3), pos(2)+pos(4)-0.25, 0.03, 0.025],'String','invert','Value',0,'BackgroundColor',get(gcf,'Color'),'Callback',@PlotCallback);
hb(2) = uicontrol('Style','checkbox','Units','normalized','Position',[pos(1)+pos(3)+0.035, pos(2)+pos(4)-0.25, 0.03, 0.025],'String','CAR','Value',0,'BackgroundColor',get(gcf,'Color'),'Callback',@PlotCallback);
N = size(DATa,1); % number of samples (minus header)
T = N/fs; % sec
set(hs,'Min',1,'Max',N,'SliderStep',[0.5/T, 2/T],'Value',1,'Callback',@PlotCallback);
gdat = guidata(gcf);
gdat.datA = DATa; gdat.datB = DATb; gdat.events = events;
gdat.fs = fs; gdat.T = T;
gdat.hs = hs; gdat.hc = hc;
gdat.he = he; gdat.hb = hb;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(~,~)
% load data
gdat = guidata(gcf);
datA = gdat.datA; datB = gdat.datB; events = gdat.events;
fs = gdat.fs; T = gdat.T;
hs = gdat.hs; hc = gdat.hc;
he = gdat.he; hb = gdat.hb;
val = round(get(hs,'Value')); % start sample
dxs = get(hc,'String');
dx = str2double(dxs{get(hc,'Value')});
ds = fix(fs*dx);
cind = val:val+ds-1;
N = length(datA);
if cind(end)>N
    cind = N-ds+1:N;
end
cdatA = datA(cind,:);
if get(hb(1),'Value') % invert A
    cdatA = -cdatA;
end
if ~isempty(datB), cdatB = datB(cind,:); end
scaleA = str2double(get(he(1),'String'));
scaleB = str2double(get(he(2),'String'));
if get(hb(2),'Value') % common average reference A and B
    cdatA = cdatA - mean(cdatA,2)*ones(1,size(cdatA,2));
    if ~isempty(datB), cdatB = cdatB - mean(cdatB,2)*ones(1,size(cdatB,2)); end
end

% events
cev = [];
if ~isempty(events)
    if isstruct(events)
        cev = events.samp(events.samp(:,1)>=cind(1) & events.samp(:,1)<=cind(end),:);
        if ~isempty(cev), cev(:,1) = cev(:,1)/fs; end
    elseif iscell(events)
        evsamp = cellfun(@(x) x*fs, events(:,2));
        cev = events(evsamp>=cind(1) & evsamp<=cind(end),:);
    end
end

% plot
cla;
t = cind/fs;
for ich = 1:size(datA,2)
    chdat = cdatA(:,ich);
    chdat = (chdat-median(chdat))/scaleA;
    plot(t,chdat+ich,'k');
    if ~isempty(datB)
        chdat = cdatB(:,ich);
        chdat = (chdat-median(chdat))/scaleB;
        plot(t,chdat+ich,'b');
    end
end
set(gca,'XLim',[val, val+ds-1]/fs);
set(hs,'SliderStep',[0.5*dx/T, 2*dx/T]);
if ~isempty(cev)
    if ~iscell(cev)
        line([cev(:,1)'; cev(:,1)'],get(gca,'YLim')'*ones(1,size(cev,1)),'Color','r');
        ctxt = events.name(cev(:,2));
        text(cev(:,1),max(get(gca,'YLim'))*ones(size(cev,1),1),ctxt,'VerticalAlignment','bottom','HorizontalAlignment','center');
    else
        for ii = 1:size(cev,1)
            line(cev{ii,2}*ones(1,2),get(gca,'YLim'),'Color','r');
            text(cev{ii,2},max(get(gca,'YLim')),cev{ii,1},'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
    end
end
end
