% cd('/Users/tnl/matlab/stereoEEG')
% parametersstats = grpstats(cc(:, 1:4), cc(:, 5), {'mean', 'std'}); 
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
nperm = 1000; shuffdims = [1 2 3 4]; % [1-Channel 2-Phase 3-Amplitude 4-Time]
[dist, permThreshold] = runPermTest(pData, shuffdims, nperm, conn);

% run test data connectivity
cc = connectivity(pData, conn, permThreshold);

%% Connected set characterization
[min, max, mean, std] = grpstats(cc(:, 1:4), cc(:, 5), {'min','max','mean','std'});

% Calculate the number of unique channels in each connected set
uniqueChannelsPerSet = arrayfun(@(set) numel(unique(connectedSets(connectedSets(:, 1) == set, 2))), unique(connectedSets(:, 1)));

% Calculate the dominant channel for each connected set
dominantChannels = arrayfun(@(set) mode(connectedSets(connectedSets(:, 1) == set, 2)), unique(connectedSets(:, 1)));

% Display results
disp('Number of Unique Channels per Connected Set:');
disp(uniqueChannelsPerSet');
disp('Dominant Channel for Each Connected Set:');
disp(dominantChannels');

%% Perform correlation analysis between phase and amplitude

uniqueSets = unique(cc(:, 5));
correlationMatrices = zeros(size(uniqueSets));

for i = 1:numel(uniqueSets)
    currentSetData = cc(cc(:, 5) == uniqueSets(i), :);
    dataForCorrelation = currentSetData(:, 2:3);
    correlationMatrix = corr(dataForCorrelation);
    correlationMatrices(i) = correlationMatrix(2);
end
%% Temporal participation ratio
figure
for s = uniqueSets'
    soi = cc(cc(:,5) == s,1:4);
    times = unique(soi(:,4));
    tcounts = zeros(size(times));
    for t = 1:length(times)
        tcounts(t) = sum(soi(:,4) == t)/length(soi);
    end
    line(times,tcounts)
    hold on
end
%% Create histogram of permutation test connectivity distribution
figure
histogram(dist, 'BinMethod','fd' ,'Normalization', 'probability');
title('Distribution of Permutation Test Connected Voxels');
xlabel('Number of Connected Voxels');
ylabel('Probability');
set(gca,'XScale','log','YScale','log')
xline(permThreshold, 'r', 'LineWidth', 1, 'Label', ['Threshold:' string(permThreshold)]);
%% Create histogram of test data connectivity distribution
figure
tdist = histcounts(cc(:,5), [unique(cc(:,5)); max(unique(cc(:,5)))+1]);
histogram(tdist, 'Normalization', 'probability','BinEdges',0:1:max(tdist) );
title('Distribution of Test Data Connected Voxels');
xlabel('Number of Connected Voxels');
ylabel('Probability');
set(gca,'XScale','log','YScale','log')
axis tight
%% Interactive 3D scatter plot of connected sets
fig = figure('Position', [100, 100, 800, 600]);
currentSet = 1;  
scatterObj = scatter3(...
    cc(cc(:,5) == currentSet, 2),...
    cc(cc(:,5) == currentSet, 3),...
    cc(cc(:,5) == currentSet, 4),...
    50, cc(cc(:,5) == currentSet, 1), 'filled');
title(['Phase Amplitude Coupling - Connected Set ' num2str(currentSet)]);
xlabel('Phase'); ylabel('Amplitude'); zlabel('Time');
[~, NA, NP, Nt] = size(pData);
xlim([0,NP+1]); ylim([0,NA+1]); zlim([0,Nt+1])
colormap(flag); colorbar('Ticks',[-1, 0, 1]);
% Plot rectangular plane at LOC
zPosition = find(PACparam.t > tones.LOC,1);
vertices = [0,0,zPosition; NP+1,0,zPosition; NP+1,NA+1,zPosition; 0,NA+1,zPosition];
faces = [1, 2, 3, 4];
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue', 'FaceAlpha', 0.25);
setappdata(fig, 'scatterObj', scatterObj);
setappdata(fig, 'sets', cc);
uniqueSets = unique(cc(:, 5));
numSets = length(uniqueSets);
setNames = cell(1, numSets);
for i = 1:numSets
    setNames{i} = ['Set ' num2str(uniqueSets(i))];
end
channelMenu = uicontrol('Style', 'popupmenu', 'Position', [10, 10, 120, 25],...
    'String', ['Select Set|' setNames], 'Callback', @updateSet);

function updateSet(source, ~)
    fig = gcf;
    scatterObj = getappdata(fig, 'scatterObj');
    cc = getappdata(fig, 'sets');

    selectedSet = source.Value - 1;  % Adjust to match the channel numbering
    if selectedSet == 0
        set(scatterObj, 'XData', cc(cc(:,5) == currentSet, 2));
        set(scatterObj, 'YData', cc(cc(:,5) == currentSet, 3));
        set(scatterObj, 'ZData', cc(cc(:,5) == currentSet, 4));
        set(scatterObj, 'CData', cc(cc(:,5) == currentSet, 1));
        title(['Interactive 3D Scatter Plot - Set ' num2str(currentSet)]);
    else
        set(scatterObj, 'XData', cc(cc(:,5) == selectedSet, 2));
        set(scatterObj, 'YData', cc(cc(:,5) == selectedSet, 3));
        set(scatterObj, 'ZData', cc(cc(:,5) == selectedSet, 4));
        set(scatterObj, 'CData', cc(cc(:,5) == selectedSet, 1));
        title(['Interactive 3D Scatter Plot - Set ' num2str(selectedSet)]);
    end
guidata(gcf);

end
%% Functions
function cc = connectivity(binMat,conn,permThreshold)
% find test data connectivity - outputs matrix of connected voxel. Each row
% corresponds to a voxel. Columns are as follows [channels, phases,
% amplitudes, times, cluster identifier]

% Find connected components in patient binary data
CC = bwconncomp(binMat,conn);
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
[channels, phases, amps, times] = ind2sub(size(binMat), pacCC);
cc = [channels, phases, amps, times, idxs];
end


function binMat = createBinaryMatrix(PACparam, PACstat, channel_labels, threshold)
% create binary matrix - outputs the binary matrix binMat containing zeros
% in all dummy channels

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


function [connDistribution, permThreshold] = runPermTest(binMat, shuffdims, nPerm, conn)
% permutation test - outputs connectivity distribution and connected set
% size threshold

cidx = find(any(binMat ~= 0, [2 3 4])); % indices of data channels 
ndim = size(binMat);
cDist = cell(nPerm,1);
for p = 1:nPerm
    permutedData = binMat;
    for n = shuffdims
        if n == 1
            idx = cidx; 
            permutedData(idx, :, :, :) = permutedData(idx(randperm(length(idx))), :, :, :);
            gctestidx = find(any(permutedData ~= 0, [2 3 4]));
            if ~isequal(gctestidx, idx)
                error('Error shuffling channel data')
            end
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