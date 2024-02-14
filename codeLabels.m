% function codeLabels()
% Generates common labeling scheme for all channel labels within seeg
% dataset
% AG 2/24

root = '/Users/tnl/Library/CloudStorage/Box-Box/sEEG';
datfile = dir(fullfile(root,'**','*_Induction.mat'));
pacchannels = cell(length(datfile),1);
nch = zeros(length(datfile),1);
for ifile = 1:length(datfile)
    load(fullfile(datfile(ifile).folder,datfile(ifile).name),'PACparam');
    disp(datfile(ifile).name);
    nch(ifile) = size(PACparam.channel_labels,1);
    pacchannels{ifile} = PACparam.channel_labels(:,2);
end

ulabels = '';
for ip = 1:length(pacchannels)
    ulabels = [ulabels; unique(pacchannels{ip})];

end
uniquelabels = unique(ulabels);
%%
% List of substrings to remove
substringsToRemove = {'right', 'left', 'ctx', 'lh', 'rh', 'cortex', 'gyrus', '-', '\s+'};

% Remove specified substrings and all whitespaces
normalized_labels = cellfun(@(s) removeSpecifiedSubstrings(lower(s), substringsToRemove), uniquelabels, 'UniformOutput', false);

function outputStr = removeSpecifiedSubstrings(inputStr, substringsToRemove)
    % Remove each specified substring
    for i = 1:length(substringsToRemove)
        inputStr = regexprep(inputStr, substringsToRemove{i}, ' ', 'ignorecase');
    end
    outputStr = inputStr;
    % % Remove all whitespace
    % outputStr = regexprep(inputStr, '\s+', '');
end
