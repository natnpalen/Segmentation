function [label_map, C] = kmeans_robust(HU_volume, mask, k, rngSeed)
%KMEANS_ROBUST A robust wrapper for k-means clustering on masked data.
%
%   This function performs k-means clustering on the finite values within a
%   masked HU volume. It includes a critical pre-check to ensure that the
%   number of valid data points is greater than the requested number of
%   clusters (k), preventing a common MATLAB error.

% --- Step 1: Prepare Data ---
% Extract the data points from within the mask
data_points = HU_volume(mask);
% Crucially, filter out any non-finite values (e.g., NaNs from interpolation)
data_points = data_points(isfinite(data_points));

% --- Step 2: Pre-flight Safety Check ---
% This is the key step. Check if we have enough data to cluster.
if numel(data_points) < k
    error('PartitionKmeans:NotEnoughData', ...
        'Cannot create %d clusters because only %d valid data points were found in the mask.', ...
        k, numel(data_points));
end

% --- Step 3: Perform K-Means Clustering ---

% --- CORRECTED: Robustly handle the random number generator seed ---
% The provided rngSeed might be a number, a string ('shuffle'), or a struct.
% The rng() function is strict. We will check if the seed is a valid number
% and fall back to the default predictable state if it is not.
if isnumeric(rngSeed) && isscalar(rngSeed) && rngSeed >= 0
    rng(rngSeed);
else
    % Fall back to the default seed if the input is not a valid number.
    % This ensures reproducible behavior.
    rng('default');
end

% Perform k-means on the 1-D vector of valid HU values
[idx, C] = kmeans(data_points, k, 'Start','plus', 'Replicates',5, 'MaxIter',1000);

% --- Step 4: Reconstruct the Label Map ---
% Create an output volume of zeros
label_map = zeros(size(HU_volume), 'uint8');
% Find the linear indices of all voxels that are both inside the mask
% AND have a finite HU value. This is where we'll place our results.
valid_indices = find(mask & isfinite(HU_volume));
% Place the cluster indices (idx) back into the correct spatial locations
label_map(valid_indices) = idx;

end