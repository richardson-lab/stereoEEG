% % parameters cd('/Users/tnl/matlab/stereoEEG')
root = '/Users/tnl/Library/CloudStorage/Box-Box/sEEG';
subj = 'HUP246_RID893'; % session name
threshold = 0.95; % for creating binary matrix
conn = conndef(4, 'minimal'); % connectivity for bwconncomp

% load data
cd(fullfile(root,subj));
load([subj '_Induction.mat']);

% create binary data matrix
pData = createBinaryMatrix(PACparam, PACstat, channel_labels, threshold);

% run permutation test
nperm = 1000; shuffdims = [1]; % [1-Channel 2-Phase 3-Amplitude 4-Time]
[dist,permThreshold] = runPermTest(pData, shuffdims, nperm, conn);

% Find connected components in patient binary data
CC = bwconncomp(pData,conn);

% Filter for connected voxel groups that are larger than perm threshold
ccN = cellfun(@numel,CC.PixelIdxList);
ccIdx = ccN > permThreshold;
pacCC = cat(1, CC.PixelIdxList{1,ccIdx});

% Create connected group identifier list
it = 1; idxs = [];
for i = ccN(ccIdx)
    ids = ones(i,1)*it;
    it = it+1;
    idxs = [idxs; ids];
end

[channels, phases, amps, times] = ind2sub(size(pData), pacCC);
ccMatrix = [channels, phases, amps, times, idxs];

%% Create histogram of permutation test connectivity distribution
figure
histogram(dist, 'BinEdges', 0:1:ceil(max(dist)), 'Normalization', 'probability');
title('Distribution of Permuatation Test Connected Voxels');
xlabel('Number of Connected Voxels');
ylabel('Probability');
xline(permThreshold, 'r', 'LineWidth', 1, 'Label', ['Threshold:' string(permThreshold)]);
%% Create histogram of patient connectivity distribution
figure
patdist = cellfun(@numel,CC.PixelIdxList);
histogram(patdist, 'Normalization', 'probability', 'BinMethod','fd');
title('Distribution of Patient Connected Voxels');
xlabel('Number of Connected Voxels');
ylabel('Probability');
set(gca,'XScale','log','YScale','log')
axis tight
%% Interactive 3D scatter plot of connected sets
fig = figure('Position', [100, 100, 800, 600]);

% filter data as needed - modify plotting code to work with filtered matrix
% Find rows with more than one unique value in column 1 for the same value in column 5
% uniqueValuesColumn5 = unique(ccMatrix(:, 5));
% filteredRows = [];
% for i = 1:length(uniqueValuesColumn5)
%     currentRows = ccMatrix(ccMatrix(:, 5) == uniqueValuesColumn5(i), :);
%     uniqueValuesColumn1 = unique(currentRows(:, 1));
%     if numel(uniqueValuesColumn1) > 1
%         filteredRows = [filteredRows; currentRows];
%     end
% end
filteredRows = ccMatrix;
% Scatter plot of significant voxels colored by channel for the first set
currentSet = 1;  
scatterObj = scatter3(...
    filteredRows(filteredRows(:,5) == currentSet, 2),...
    filteredRows(filteredRows(:,5) == currentSet, 3),...
    filteredRows(filteredRows(:,5) == currentSet, 4),...
    50, filteredRows(filteredRows(:,5) == currentSet, 1), 'filled');
title(['Interactive 3D Scatter Plot - Set ' num2str(currentSet)]);
xlabel('Phase'); ylabel('Amplitude'); zlabel('Time');
[~, NA, NP, Nt] = size(pData);
xlim([0,NP+1]); ylim([0,NA+1]); zlim([0,Nt+1])
colormap(flag); colorbar;

% Store scatterObj and sets as appdata of the figure
setappdata(fig, 'scatterObj', scatterObj);
setappdata(fig, 'sets', filteredRows);

% Create a uicontrol dropdown menu for channel selection
uniqueSets = unique(filteredRows(:, 5));
numSets = length(uniqueSets);
setNames = cell(1, numSets);
for i = 1:numSets
    setNames{i} = ['Set ' num2str(uniqueSets(i))];
end

channelMenu = uicontrol('Style', 'popupmenu', 'Position', [10, 10, 120, 25],...
    'String', ['Select Set|' setNames], 'Callback', @updateSet);

% Function to update the scatter plot based on the selected channel
function updateSet(source, ~)
    fig = gcf;
    scatterObj = getappdata(fig, 'scatterObj');
    filteredRows = getappdata(fig, 'sets');

    selectedSet = source.Value - 1;  % Adjust to match the channel numbering
    if selectedSet == 0
        % Show default channel
        set(scatterObj, 'XData', filteredRows(filteredRows(:,5) == currentSet, 2));
        set(scatterObj, 'YData', filteredRows(filteredRows(:,5) == currentSet, 3));
        set(scatterObj, 'ZData', filteredRows(filteredRows(:,5) == currentSet, 4));
        set(scatterObj, 'CData', filteredRows(filteredRows(:,5) == currentSet, 1));
        title(['Interactive 3D Scatter Plot - Set ' num2str(currentSet)]);
    else
        % Show the selected channel
        set(scatterObj, 'XData', filteredRows(filteredRows(:,5) == selectedSet, 2));
        set(scatterObj, 'YData', filteredRows(filteredRows(:,5) == selectedSet, 3));
        set(scatterObj, 'ZData', filteredRows(filteredRows(:,5) == selectedSet, 4));
        set(scatterObj, 'CData', filteredRows(filteredRows(:,5) == selectedSet, 1));
        title(['Interactive 3D Scatter Plot - Set ' num2str(selectedSet)]);
    end
guidata(gcf);

end
%%

% create binary matrix - outputs the binary matrix binMat containing zeros
% in all dummy channels
function binMat = createBinaryMatrix(PACparam, PACstat, channel_labels, threshold)

gcidx = find(ismember(channel_labels(:,1),PACparam.channel_labels(:,1)));
bcidx = find(~ismember(channel_labels(:,1),PACparam.channel_labels(:,1)));
Tch = length(channel_labels);
[~, NP, NA, Nt, ~] = size(PACstat);
binMat = zeros(Tch, NP, NA, Nt);

for ich = 1:Tch
    if any(ismember(gcidx,ich))
        binMat(ich,:,:,:) = PACstat(gcidx == ich,:,:,:,1) >= threshold;
    elseif any(ismember(bcidx,ich))
        binMat(ich,:,:,:) = zeros(NP,NA,Nt);
    end
end
end

% permutation test - outputs connectivity distribution and connected set
% size threshold
function [connDistribution, permThreshold] = runPermTest(binMat, shuffdims, nPerm, conn)
cidx = find(any(binMat ~= 0, [2 3 4])); % indices of data channels 
ndim = size(binMat);
cDist = cell(nPerm,1);
for p = 1:nPerm
    permutedData = binMat;
    for n = shuffdims
        if n == 1
            idx = cidx; 
            permutedData(idx, :, :, :) = permutedData(idx(randperm(length(idx))), :, :, :);
            % gctestidx = find(any(permutedData ~= 0, [2 3 4]));
            % isequal(gctestidx, gcidx);
        elseif n == 2
            idx = 1:ndim(n); 
            permutedData(:, idx, :, :) = permutedData(:, idx(randperm(length(idx))), :, :);
        elseif n == 3
            idx = 1:ndim(n); 
            permutedData(:, :, idx, :) = permutedData(:, :, idx(randperm(length(idx))), :);
        elseif n == 4
            idx = 1:ndim(n); 
            permutedData(:, :, :, idx) = permutedData(:, :, :, idx(randperm(length(idx))));
        end
    end
    CC = bwconncomp(permutedData,conn); % determine connectedness
    cDist{p} = CC.PixelIdxList; % store connected voxel matrix
end
connDistributionCell = cellfun(@(nestedCellArray) cellfun(@numel, nestedCellArray), cDist, 'UniformOutput', false);
connDistribution = horzcat(connDistributionCell{:});
permThreshold = quantile(connDistribution,0.95);
end