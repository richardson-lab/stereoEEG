function seegCFbatch
% function seegCFbatch
%   Distribution of eigenvalue criticiality x frequency obtained from
%   stability analysis of sEEG data using multivariate autoregression
%   modeling: batch analysis (looping through all available sessions).
%
%   DR 06/2022

% parameters
ddir = '/Users/drew/Box/sEEG/'; % data directory
dprint = '/Users/drew/Desktop/'; % directory to print figures
tstate = {'HUP223_Emergence',[25 32],[0 7]; % subject, [start stop] awake, [start stop] anesthetized (in minutes)
          'HUP224_Emergence',[15 22],[0 7];
          'HUP224_Induction',[0 10],[20 30];
          'HUP225_Emergence',[15 20],[0 5];
          'HUP225_Induction',[0 10],[15 25];
          'HUP227_Emergence',[8 13],[0 5];
          'HUP227_Induction',[0 10],[30 40];
          'HUP229_Emergence',[7 12],[0 5]};
      
% main
cd(ddir);
for ii = 1:size(tstate,1)
    subj = tstate{ii,1};
    tawake = tstate{ii,2};
    tanesth = tstate{ii,3};
    load(fullfile(subj(1:6),subj,[subj '.mat']),'session');
    disp(subj);
    seegCF(session,tawake,tanesth);
    drawnow;
    orient landscape
    print(fullfile(dprint,[mfilename '.ps']),'-dpsc2','-fillpage','-append');
    clear session
end
