%-------------------------------------------------------------------------------------
%plot the figures
load('resultsindividual');
 t=1:size(p,2)-1; 

 figure;

 
 %plot learning curve
 subplot(211);  
 h=plot(t, pmid(2:end),'b-'); set(h, 'LineWidth',2);
 hold on;
 plot(t, p05(2:end),'b', t, p95(2:end), 'b');
 if(MaxResponse == 1)
      hold on; [y, x] = find(Responses > 0);
      h = plot(x,y,'o'); set(h, 'MarkerFaceColor',[.9 .9 .9]);
      set(h, 'MarkerEdgeColor', 'k' ,'MarkerSize', 4);
      hold on; [y, x] = find(Responses == 0);
      h = plot(x,zeros(size(x)),'o'); set(h, 'MarkerFaceColor', [.9 .9 .9]);
      set(h, 'MarkerEdgeColor', 'k','MarkerSize', 4);
      axis([1 t(end)  0 1.0]);
 else
      hold on; plot(t, Responses./MaxResponse,'ko');
      axis([1 t(end)  0 1]);
 end
 line([1 t(end)], [BackgroundProb  BackgroundProb ]);
 title(['IO(0.95) Learning trial = ' num2str(cback) '   Learning state process variance = ' num2str(SigE^2) ]);
 xlabel('Trial Number')
 ylabel('Probability of a Correct Response')

%  %plot IO certainty
 subplot(223)
 plot(t,pcert(2:end),'k')
 line([ 1 t(end)],[0.90 0.90]);
 line([ 1 t(end)],[0.99 0.99]);
 line([ 1 t(end)],[0.95 0.95]);
 axis([1 t(end)  0 1]);
 grid on;
 xlabel('Trial Number')
 ylabel('Certainty')
 