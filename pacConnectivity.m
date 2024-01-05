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
shuffdims = [2 3]; % [1-Channel 2-Phase 3-Amplitude 4-Time]
nperm = 1000;
[dist,permThreshold] = runPermTest(pData, shuffdims, nperm, conn);

% Find connected components in patient binary data
CC = bwconncomp(pData,conn);

% Filter for connected voxel groups that are larger than perm threshold
ccIdx = cellfun(@numel, CC.PixelIdxList) > permThreshold;
pacCC = cat(1, CC.PixelIdxList{1,ccIdx});
[sizeC, sizeP, sizeA, sizeT] = size(pData);
[channels, phases, amps, times] = ind2sub([sizeC, sizeP, sizeA, sizeT], pacCC);
rowSubscriptsMatrix = [channels, phases, amps, times];

%% Create a 3D scatter plot
figure
scatter3(rowSubscriptsMatrix(:, 2), rowSubscriptsMatrix(:, 3), rowSubscriptsMatrix(:, 4), 50, rowSubscriptsMatrix(:, 1), 'filled');
title([subj(1:6) ' Connected Voxels']);
xlabel('Phase');
ylabel('Amplitude');
zlabel('Time');
colormap(jet);
cbar = colorbar;
ticks = unique(rowSubscriptsMatrix(:, 1));
cbar.Ticks = ticks;
cbar.TickLabels = int2str(ticks);

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