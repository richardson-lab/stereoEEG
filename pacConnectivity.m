% % parameters
root = '/Users/andrewgabros/Library/CloudStorage/Box-Box/sEEG';
subj = 'HUP246_RID893'; % session name
threshold = 0.95;
% load data
cd(fullfile(root,subj));
load([subj '_Induction.mat']);

% define connectivity for analysis
conn = conndef(4, 'minimal');

% create binary data matrix
[pData, gcidx] = createBinaryMatrix(PACparam, PACstat, channel_labels, threshold);

% run permutation test
[dist,permThreshold] = runPermTest(pData, gcidx, 1000, conn);

% Find connected components in patient binary data
CC = bwconncomp(pData,conn);

% Filter for connected voxel groups that are larger than perm threshold
ccIdx = cellfun(@numel, CC.PixelIdxList) > permThreshold;
% pacCC = CC.PixelIdxList{(cellfun(@numel, CC.PixelIdxList) > permThreshold)};
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
histogram(patdist, 'Normalization', 'probability');
title('Distribution of Patient Connected Voxels');
xlabel('Number of Connected Voxels');
ylabel('Probability');
% xline(permThreshold, 'r', 'LineWidth', 1, 'Label', ['Threshold:' string(permThreshold)]);
set(gca,'XScale','log','YScale','log')
axis tight
%%
function [binMat, gcidx] = createBinaryMatrix(PACparam, PACstat, channel_labels, threshold)

gcidx = find(ismember(channel_labels(:,1),PACparam.channel_labels(:,1)));
bcidx = find(~ismember(channel_labels(:,1),PACparam.channel_labels(:,1)));
Tch = length(channel_labels);
[~, NP, NA, Nt, ~] = size(PACstat);
binMat = zeros(Tch, NP, NA, Nt);

for ich = 1:Tch
    if any(ismember(gcidx,ich))
        binMat(ich, :,:,:) = PACstat(gcidx == ich,  :,:,:,1) >= threshold;
    elseif any(ismember(bcidx,ich))
        binMat(ich,  :,:,:) = zeros(NP,NA,Nt);
    end
end
end

% permutation test - need to add capability to shuffle other dimensions

function [connDistribution, permThreshold] = runPermTest(binMat, gcidx, nPerm, conn)

[Nch, NP, NA, Nt, ~] = size(binMat);
cDist = cell(nPerm,1);
for p = 1:nPerm
    % Reshape the data matrix to a 2D matrix where each column corresponds to a channel
    reshapedData = reshape(binMat, [], Nch);

    % Randomly permute only the specified channels
    permutedData = zeros(size(reshapedData));
    permutedData(:,gcidx) = reshapedData(:, gcidx(randperm(length(gcidx))));

    % Reshape the permuted data back to the original dimensions
    permutedData = reshape(permutedData, Nch, NA, NP, Nt);

    % determine connectedness
    CC = bwconncomp(permutedData,conn);

    % store connected voxel matrix
    cDist{p} = CC.PixelIdxList;
end
connDistributionCell = cellfun(@(nestedCellArray) cellfun(@numel, nestedCellArray), cDist, 'UniformOutput', false);
connDistribution = horzcat(connDistributionCell{:});
permThreshold = quantile(connDistribution,0.95);
end