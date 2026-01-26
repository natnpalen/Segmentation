function [labels_raw, centroids_sorted] = kmeans_label_masked(HU, mask, k, rngSeed)
% KMEANS_LABEL_MASKED  K-means on in-mask HU; relabel 1..k by mean HU (low→high).
% Returns labels volume (uint8), 0 outside mask.

if nargin<4, rngSeed = []; end
vals = double(HU(mask));
if isempty(vals)
    labels_raw = zeros(size(HU), 'uint8');
    centroids_sorted = zeros(k,1);
    return;
end

% Optional reproducibility
if ~isempty(rngSeed), rng(rngSeed); end

% Run K-means
idx = kmeans(vals(:), k, 'Replicates', 3, 'MaxIter', 200);

% Compute centroids and order low→high
centroids = zeros(k,1);
for ii = 1:k
    centroids(ii) = mean(vals(idx==ii));
end
[~, order] = sort(centroids, 'ascend');

% Relabel to 1..k by ordered centroids
newL = zeros(size(idx));
for ii = 1:k
    newL(idx == order(ii)) = ii;
end

% Rasterize
labels_raw = zeros(size(HU), 'uint8');
labels_raw(mask) = uint8(newL);
centroids_sorted = sort(centroids, 'ascend');
end
