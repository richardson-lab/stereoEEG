function seegStabilityBatch
% function seegStabilityBatch
%   Stability analysis of sEEG data using multivariate autoregression
%   modeling: batch analysis (looping through all available sessions).
%
%   DR 03/2022

% parameters
ddir = '/Users/drew/Box/sEEG/'; % data directory

% main
cd(ddir);
fsubj = dir('HUP*');
Nsubj = length(fsubj);
for ii = 1:Nsubj
    try % emergence session
        load(fullfile(fsubj(ii).name,[fsubj(ii).name '_Emergence'],[fsubj(ii).name '_Emergence.mat']),'session');
        disp([fsubj(ii).name '_Emergence']);
        seegStability(session);
        drawnow;
        clear session
    end
    try % induction session
        load(fullfile(fsubj(ii).name,[fsubj(ii).name '_Induction'],[fsubj(ii).name '_Induction.mat']),'session');
        disp([fsubj(ii).name '_Induction']);
        seegStability(session);
        drawnow;
        clear session
    end
end
