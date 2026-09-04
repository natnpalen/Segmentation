function opts = quality_defaults(opts)
% QUALITY_DEFAULTS  Thresholds used to judge a bone mask and to decide
% whether a fallback segmentation pass is worth running.
%
%   opts = bone.quality_defaults(opts)
%
% Fills in any of these that are missing, leaving values the caller set:
%
%   LowDensityHU    : 250   interior median below this is a low-density
%                           (osteoporotic) bone — flagged, not rejected
%   MinContrastHU   : 150   interior median must sit at least this far above
%                           the surrounding non-air material, otherwise the
%                           mask boundary is not backed by a density step
%   MaxLowFrac      : 0.35  fraction of mask voxels allowed below
%                           TissueCeilingHU before tissue is suspected
%   TissueCeilingHU : 150   HU below which a voxel is not clearly bone
%   MinFillRatio    : 0.50  mask volume as a fraction of the blob it grew
%                           from; below this the bone is under-segmented
%   SourceVolMM3    : 0     volume of that source blob (0 = unknown)

DEFAULTS = struct( ...
    'LowDensityHU',    250, ...
    'MinContrastHU',   150, ...
    'MaxLowFrac',      0.35, ...
    'TissueCeilingHU', 150, ...
    'MinFillRatio',    0.50, ...
    'SourceVolMM3',    0 ...
);

names = fieldnames(DEFAULTS);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = DEFAULTS.(names{i});
    end
end
end
