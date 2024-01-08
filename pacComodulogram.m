function pacComodulogram(tones,PACparam,MI,MIpval)
% function pacComodulogram(tones,PACparam,MI,MIpval)
%   Plot PAC comodulogram + behavior vs time for sEEG data. Inputs come
%   from varaibles saved to .mat file using pacMI.m and toneResponseAnaly.m
%
%   DR 05/2023

% behavior vs time
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/12 1/6 2/3 2/3],'Color','w');
ax1 = subplot(3,4,1:3);
patch([tones.ttrial', fliplr(tones.ttrial')],[tones.p95 fliplr(tones.p05)],'k','FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
hold on; plot(tones.ttrial,tones.pmid,'k','LineWidth',3);
ind = find(tones.responses);
plot(tones.ttrial(ind),ones(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
ind = find(~tones.responses);
plot(tones.ttrial(ind),zeros(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[tones.ttrial(1) tones.ttrial(end)],'YLim',[0 1]);
line(get(gca,'XLim'),0.05*ones(1,2),'Color','k','LineStyle','--','LineWidth',2);
line(tones.LOC*ones(1,2),get(gca,'YLim'),'Color','r','LineStyle','-','LineWidth',2);
text(tones.LOC,max(get(gca,'YLim')),[' LOC = ' num2str(round(tones.LOC)) ' s'],'HorizontalAlignment','left','VerticalAlignment','top','Color','r','FontSize',14);
xlabel('time (s)');
ylabel('response probability');

% threshold MI based on significance
if ~isempty(MIpval)
    ind = find(MIpval > .05);
    MI(ind) = 0;
end

% comodulogram vs time: projection 1
[NP,NA,Nt] = size(MI);
colormap(flipud(hot(512)));
ax2 = subplot(3,4,5:7);
rproj = squeeze(max(MI,[],2));
imagesc(PACparam.t, 1:NP, rproj);
set(gca,'Box','off','YDir','normal','TickDir','out','FontSize',14,'CLim',[0 .02]);
ytick = get(gca,'YTick'); set(gca,'YTickLabel',num2str(mean(PACparam.rP(ytick,:),2),'%4.2f\n'));
ylabel('phase freq (Hz)');
hc = colorbar('northoutside');
set(get(hc,'YLabel'),'String','phase-amplitude coupling modulation index');

% comodulogram vs time: projection 2
ax3 = subplot(3,4,9:11);
rproj = squeeze(max(MI,[],1));
imagesc(PACparam.t, 1:NA, rproj);
set(gca,'Box','off','YDir','normal','TickDir','out','FontSize',14,'CLim',[0 .02]);
ytick = get(gca,'YTick'); set(gca,'YTickLabel',num2str(mean(PACparam.rA(ytick,:),2),'%4.2f\n'));
xlabel('time (s)'); ylabel('amplitude freq (Hz)');
xlim1 = get(ax1,'XLim');
xlim2 = get(ax2,'XLim');
xlim3 = get(ax3,'XLim');
xlim = [max([xlim1(1), xlim2(1), xlim3(1)]), min([xlim1(2), xlim2(2), xlim3(2)])];
set(ax1,'XLim',xlim);
set(ax2,'XLim',xlim);
set(ax3,'Xlim',xlim);

% interactive comodulogram plot
pos = get(ax3,'Position');
[~,imin] = min(abs(PACparam.t-xlim(1)));
[~,imax] = min(abs(PACparam.t-xlim(2)));
hs = uicontrol('Style','slider','Units','normalized','Position',[pos(1), pos(2)+pos(4), pos(3), 0.025],'Min',imin,'Max',imax,'SliderStep',[1/(imax-imin), 5/(imax-imin)],'Value',imin,'Callback',@PlotCallback);
ax4 = subplot(3,4,[8 12]);
gdat = guidata(gcf);
gdat.PACparam = PACparam; gdat.MI = MI;
gdat.ax2 = ax2; gdat.ax3 = ax3; gdat.ax4 = ax4;
gdat.hs = hs;
guidata(gcf,gdat);
PlotCallback(hs);
end

function PlotCallback(~,~)
gdat = guidata(gcf);
iT = round(gdat.hs.Value);
cMI = squeeze(gdat.MI(:,:,iT));
axes(gdat.ax4); cla;
fP = mean(gdat.PACparam.rP,2);
fA = mean(gdat.PACparam.rA,2);
imagesc(1:length(fP),1:length(fA),cMI);
hc = colorbar('northoutside'); set(get(hc,'YLabel'),'String','modulation index');
axis square; set(gca,'Box','off','YDir','normal','YAxisLocation','right','TickDir','out','FontSize',14);
xtick = get(gca,'XTick'); ytick = get(gca,'YTick'); set(gca,'XTickLabel',num2str(round(fP(xtick),2),'%4.2f\n'),'YTickLabel',num2str(round(fA(ytick),1),'%5.1f\n'));
xlabel('phase frequency (Hz)'); ylabel('amplitude frequency (Hz)');
axes(gdat.ax2);
delete(findobj(gdat.ax2,'Type','line'));
line(gdat.PACparam.t(iT)*ones(1,2),get(gca,'YLim'),'Color','k','LineWidth',1);
axes(gdat.ax3);
delete(findobj(gdat.ax3,'Type','line'));
line(gdat.PACparam.t(iT)*ones(1,2),get(gca,'YLim'),'Color','k','LineWidth',1);
end
