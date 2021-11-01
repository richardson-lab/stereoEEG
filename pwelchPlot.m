function [Pall,f] = pwelchPlot(dat,fs,varargin)
% function pwelchPlot(dat,fs,frq,L,O,N,chn,hax)
%   Plot of Welch's power spectral density estimate.
%
%   dat = data vector (samples x channels)
%   fs = sampling frequency (Hz)
%   frq = frequncy range [min max] (Hz)
%   L = window length (s)
%   O = window overlap (s)
%   N = number of logorithmically-spaced power estimates over frequency range
%   chn = channel numbers/names
%   hax = axis handle for plot
% 
%   DR 03/2018

% parameters
frq = [0.4 300]; % frequency range (Hz)
L = 5; % window length (s)
O = 1; % window overlap (s)
N = 500; % number of logorithmically-spaced power estimates over frequency range
chn = []; % channel numbers/names
if nargin>=3, frq = varargin{1}; end
if nargin>=4, L = varargin{2}; end
if nargin>=5, O = varargin{3}; end
if nargin>=6, N = varargin{4}; end
if nargin>=7, chn = varargin{5}; end
if nargin>=8
    hax = varargin{6};
    axes(hax);
else
    figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
end

% analysis
Nch = size(dat,2); % channels
T = round(size(dat,1)/fs); % s
L = L*fix(fs); % samples
O = O*fix(fs); % samples
if Nch>1, clr = colormap(lines(Nch)); else, clr = [0 0 0]; end
for ich = 1:Nch
    cdat = dat(:,ich);
    cdat = cdat-median(cdat);
    [P,f] = pwelch(cdat,L,O,logspace(log10(frq(1)),log10(frq(2)),N),fs);
    ind = find(f>=frq(1) & f<=frq(2));
    f = f(ind); P = P(ind);
    if ~isempty(chn) && length(chn)==Nch
        if iscell(chn)
            plot(f,P,'Color',clr(ich,:),'Tag',chn{ich}); hold on;
        else
            plot(f,P,'Color',clr(ich,:),'Tag',num2str(chn(ich))); hold on;
        end
    else
        plot(f,P,'Color',clr(ich,:),'Tag',num2str(ich)); hold on;
    end
    if ich==1
        Pall = zeros(Nch,length(P));
    end
    Pall(ich,:) = P;
    if ich == 1
        set(gca,'Box','off','XScale','log','YScale','log','XLim',frq,'TickDir','out','FontSize',14); axis tight; grid on;
        xlabel('frequency (Hz)'); ylabel('power (\muV^2/Hz)');
    end
    drawnow
end
text(max(get(gca,'XLim')),max(get(gca,'YLim')),[num2str(T) '-s recording'],'HorizontalAlignment','right','VerticalAlignment','top');
% axis square;
