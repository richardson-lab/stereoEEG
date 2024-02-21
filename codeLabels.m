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
    pacchannels{ifile} = PACparam.channel_labels;
end
save('ChannelLabels.mat','pacchannels');

ulabels = '';
for ip = 1:length(pacchannels)
    ulabels = [ulabels; unique(pacchannels{ip}(:,2))];
end
uniquelabels = unique(ulabels);
%%
% List of substrings to remove
substrings_to_remove = {'right', 'left', 'ctx', 'lh', 'rh', 'cortex', 'gyrus', '-'};

% Normalize labels - lowercase, remove substrings, hyphens and whitespace
normalized_labels = cellfun(@(s) normalizeLabels(lower(s), substrings_to_remove), uniquelabels, 'UniformOutput', false);

intlabels = groupSimilarStrings(normalized_labels);

map = cell(max(intlabels),1);
for i = 1:length(map)
    map{i} = uniquelabels(intlabels==i);
end

save('ChannelLabelsMap.mat','map');

%% FUNCTIONS %%

function outputStr = normalizeLabels(inputStr, substrings_to_remove)
    % Remove each specified substring
    for i = 1:length(substrings_to_remove)
        inputStr = erase(inputStr, substrings_to_remove{i});
    end

    % Remove all whitespace
    outputStr = regexprep(inputStr, '\s+', '');
end

function labels = groupSimilarStrings(strings)
    n = numel(strings);
    labels = zeros(n, 1);
    groupCount = 0;
    normalizedStrings = lower(strings);
    for i = 1:n
        if labels(i) == 0
            groupCount = groupCount + 1; % New group
            labels(i) = groupCount;
        end
        
        for j = i+1:n
            if labels(j) == 0 && areSimilar(normalizedStrings{i}, normalizedStrings{j})
                labels(j) = labels(i); % Assign to the same group
            end
        end
    end
end

function isSimilar = areSimilar(str1, str2)
    % Define a simple criterion for similarity
    isSimilar = contains(str1, str2) || contains(str2, str1);
end
