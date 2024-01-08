% Comparing PAC significance based on permutation test of Mukamel et al
% 2014 J Neurosci vs. von Mises gof.
%
% DR 09/2023

% load('/Users/drew/Library/CloudStorage/Box-Box/sEEG/HUP225_RID700/HUP225_RID700_Induction.mat')
% mi1 = PACmi.RC01;
% pval1 = PACmi.RC01pval;
% ich = find(matches(PACparam.channel_labels(:,1),'RC01'));
% mi2 = squeeze(PACstat(ich,:,:,:,3));
% pval2 = 1-squeeze(PACstat(ich,:,:,:,1));

% load('/Users/drew/Library/CloudStorage/Box-Box/sEEG/HUP230_RID808/HUP230_RID808_Induction.mat')
% mi1 = PACmi.LB03;
% pval1 = PACmi.LB03pval;
% ich = find(matches(PACparam.channel_labels(:,1),'LB03'));
% mi2 = squeeze(PACstat(ich,:,:,:,3));
% pval2 = 1-squeeze(PACstat(ich,:,:,:,1));

load('/Users/drew/Library/CloudStorage/Box-Box/sEEG/HUP239_RID847/HUP239_RID847_Induction.mat')
mi1 = PACmi.RB01;
pval1 = PACmi.RB01pval;
ich = find(matches(PACparam.channel_labels(:,1),'RB01'));
mi2 = squeeze(PACstat(ich,:,:,:,3));
pval2 = 1-squeeze(PACstat(ich,:,:,:,1));

% load('/Users/drew/Library/CloudStorage/Box-Box/sEEG/HUP241_RID890/HUP241_RID890_Induction.mat')
% mi1 = PACmi.RB01;
% pval1 = PACmi.RB01pval;
% % pacComodulogram(tones,PACparam,mi1,pval1);
% ich = find(matches(PACparam.channel_labels(:,1),'RB01'));
% mi2 = squeeze(PACstat(ich,:,:,:,3));
% pval2 = 1-squeeze(PACstat(ich,:,:,:,1));
% % pacComodulogram(tones,PACparam,mi2,pval2);

% compare
figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
h = histogram(log10(mi1)); hold on;
ind1 = find(pval1<.05);
histogram(log10(mi1(ind1)),h.BinEdges);
ind2 = find(pval2<.05);
histogram(log10(mi2(ind2)),h.BinEdges);
set(gca,'Box','off','TickDir','out','FontSize',12)
xlabel('log_{10}(modulation index)');
ylabel('count');
legend({'all','perm p < .05','gof p < .05'});
twoinone = length(intersect(ind1,ind2))/length(ind2)*100;
title([num2str(twoinone,'%4.1f') '% gof are same as perm']);
