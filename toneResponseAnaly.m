function toneResponseAnaly
% function toneResponseAnaly
%   Analysis of tone response data (run toneResponseSync.m first). Computes
%   time-continuous probability of response based on approach of Smith et
%   al. (2004)  Dynamic Analysis of Learning in Behavioral Experiments. The
%   Journal of Neuroscience 24(2):447-461.
%
%   DR 12/2022

% parameters
ddir = '/Users/drew/Box/sEEG/'; % root data directory
subj = 'HUP250_RID534'; % subject (leave empty to batch process all subjects)
thresh = 2; % reaction time threshold for tone response (s)
        
% load tone response data
cd(ddir);
if isempty(subj)
    fsubj = dir('HUP*');
    Nsubj = length(fsubj);
else
    fsubj(1).name = subj;
    Nsubj = 1;
end
for ii = 1:Nsubj
    try % tone response task only happend on induction sessions
        load(fullfile(fsubj(ii).name,[fsubj(ii).name '_Induction.mat']),'annotations','t');
        disp([fsubj(ii).name '_Induction']);
    catch
        continue;
    end
    
    % compute reaction time (RT) to each tone (inf if no response)
    itrial = strcmp(annotations(:,1),'tone on');
    ttrial = cell2mat(annotations(itrial,2));
    if t(end)-ttrial(end) > 15 % fill in tone on events through end of session since task often stopped after LOC clearly obtained
        tint = [5 15]; % same inter-stimulus time interval range (s) used in toneResponse.m
        ct = ttrial(end);
        while ct < t(end)
            isi = tint(1) + diff(tint)*rand;
            ttrial(end+1) = ct + isi;
            ct = ct + isi;
        end
    end
    iresp = strcmp(annotations(:,1),'button press');
    tresp = cell2mat(annotations(iresp,2));
    RT = Inf*ones(1,length(ttrial)); % reaction time
    for jj = 1:length(ttrial)
        deltaR = tresp-ttrial(jj);
        deltaR(deltaR<=0) = [];
        deltaT = ttrial-ttrial(jj);
        deltaT(deltaT<=0) = [];
        if isempty(deltaR) || (deltaT(1)<deltaR(1)) % if no following button press or next tone occurs before next response, then no response
            continue;
        end
        RT(jj) = deltaR(1); % first button press after tone
    end
    
    % binarize responses based on RT threshold
    responses = RT<thresh;
    
    % binomial expectation maximization analysis (Smith et al 2004)
    runanalysis(responses,1,0.5);
    
    % loss of consciousness (LOC) criteria (Purdon et al 2013 PNAS)
    load('resultsindividual','pmid','p05','p95');
    pmid = pmid(2:end); p05 = p05(2:end); p95 = p95(2:end);
    ind = find(pmid < .05);
    if ~isempty(ind)
        LOC = ttrial(ind(1));
    else
        LOC = ttrial(end);
    end
    
    % plot results
    figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    patch([ttrial', fliplr(ttrial')],[p95 fliplr(p05)],'k','FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
    hold on; plot(ttrial,pmid,'k','LineWidth',3);
    ind = find(responses);
    plot(ttrial(ind),ones(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
    ind = find(~responses);
    plot(ttrial(ind),zeros(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[ttrial(1) ttrial(end)],'YLim',[0 1]);
    line(get(gca,'XLim'),0.05*ones(1,2),'Color','k','LineStyle','--','LineWidth',2);
    line(LOC*ones(1,2),get(gca,'YLim'),'Color','r','LineStyle','-','LineWidth',2);
    text(LOC,max(get(gca,'YLim')),[' LOC = ' num2str(round(LOC)) ' s'],'HorizontalAlignment','left','VerticalAlignment','top','Color','r','FontSize',14);
    xlabel('time (s)')
    ylabel('response probability')
    title([fsubj(ii).name '_Induction'],'Interpreter','none');
    drawnow;
    orient landscape
    print toneResponseAnaly.ps -dpsc2 -fillpage -append
    if Nsubj > 1
        pause(5); close;
    end
end
