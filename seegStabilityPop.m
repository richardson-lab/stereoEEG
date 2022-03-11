function seegStabilityPop(stability)
% function seegStabilityPop(stability)
%   Stability analysis of sEEG data using multivariate autoregression
%   modeling: population analysis.
%
%   Plots change in fraction of modes with lambda above a threshold value
%   over time, aligned on annotated induction/emergence events. stability =
%   structured saved in seegStability.mat by seegStability.m
%
%   DR 03/2022

% parameters
tbaseline = 300; % time at beginning of session to use as baseline (s)

% prepare plot
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
h1 = subplot(1,2,1); hold on; set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('\Delta criticality (%)'); title('Emergence'); xlabel('time (min)');
h2 = subplot(1,2,2); hold on; set(gca,'Box','off','TickDir','out','FontSize',14);
ylabel('\Delta criticality (%)'); title('Induction'); xlabel('time (min)');
Nsess = length(stability);
clr = colormap(parula(Nsess));
for ii = 1:Nsess
    
    % relative changes
    dt = (stability(ii).t-tbaseline)/60; % relative time (min)
    lambthresh = stability(ii).lambthresh;
    lbaseline = mean(lambthresh(dt<=0)); % average modes above threshold in baseline period
    lambthresh = (lambthresh-lbaseline)*100; % percent change in modes above threshold
    
    % plot
    if contains(stability(ii).subj,'Emergence')
        plot(h1,dt,lambthresh,'Color',clr(ii,:),'LineWidth',3,'Tag',stability(ii).subj);
    else
        plot(h2,dt,lambthresh,'Color',clr(ii,:),'LineWidth',3,'Tag',stability(ii).subj);
    end
end
axis(h1,'tight','square');
axis(h2,'tight','square');
ylim(1,:) = get(h1,'YLim');
ylim(2,:) = get(h2,'YLim');
set(h1,'YLim',[min(ylim(:,1)), max(ylim(:,2))],'YGrid','on');
set(h2,'YLim',[min(ylim(:,1)), max(ylim(:,2))],'YGrid','on');
