function labels = project_labels_hires(HU_hires, mask_hires, centroids_sorted, hu_min, hu_max)
% PROJECT_LABELS_HIRES
% Deterministically project hi-res voxels to k density classes using
% low-res k-means centroids. Produces labels 1..k (low->high), 0 outside mask.
%
% Inputs
%   HU_hires          : 3D double (NaN outside mask is fine)
%   mask_hires        : 3D logical
%   centroids_sorted  : kx1 double, strictly ascending (low->high)
%   hu_min, hu_max    : doubles (robust clamp from low-res; optional but recommended)
%
% Output
%   labels            : uint8, 0 outside mask, 1..k inside mask

assert(isvector(centroids_sorted) && all(isfinite(centroids_sorted)), ...
    'centroids_sorted must be a finite vector');
centroids_sorted = double(centroids_sorted(:));
k = numel(centroids_sorted);
assert(k >= 2, 'Need at least 2 centroids.');

HU = double(HU_hires);
M  = logical(mask_hires);

% Clamp (optional but stabilizes tails)
if nargin >= 4 && ~isempty(hu_min) && ~isempty(hu_max) && isfinite(hu_min) && isfinite(hu_max)
    HU = max(min(HU, hu_max), hu_min);
end

% Build k-1 thresholds as midpoints between centroids
T = 0.5 * (centroids_sorted(1:end-1) + centroids_sorted(2:end));

% Edges for discretization: [-Inf, T..., +Inf]
edges = [-inf; T; +inf];

% Prepare output
labels = zeros(size(HU), 'uint8');

% Classify only valid in-mask voxels
valid = M & isfinite(HU);
vals  = HU(valid);

% Bin using discretize: returns 1..k where edges are half-open intervals
bin = discretize(vals, edges);  % 1..k

% Safety: ensure no zeros
bin(~isfinite(vals)) = 1; % shouldn't happen due to 'valid', but guard anyway

labels(valid) = uint8(bin);
labels(~M) = 0;
end
