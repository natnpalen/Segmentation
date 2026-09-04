function q = mask_quality(vol, mask, spacing, opts)
% MASK_QUALITY  Judge whether a bone mask is bone, and how confidently.
%
%   q = bone.mask_quality(ds.HU, mask, ds.spacing, opts)
%
% Measures the mask against the material immediately around it rather than
% against absolute HU numbers, so it stays meaningful for an osteoporotic
% bone whose every value is low. The question is not "is this dense?" but
% "is this denser than what surrounds it, and is what we kept consistent
% with itself?".
%
% Fields
%   core_hu     median HU of the mask interior (>1 mm from the surface)
%   surface_hu  median HU of the outermost voxel shell
%   rind_hu     median HU of non-air material just outside the mask —
%               the soft tissue we are trying not to grade
%   contrast    core_hu - rind_hu. Low contrast means the mask boundary is
%               not backed by a real density step.
%   low_frac    fraction of mask voxels below TissueCeilingHU
%   volume_mm3  mask volume
%   fill_ratio  volume_mm3 / source_vol_mm3 when the source blob volume is
%               known (0 otherwise) — well under 1 means under-segmented
%   score       0..1 summary used to compare two candidate masks
%   flags       cellstr, any of 'low_density', 'tissue_suspect',
%               'under_segmented', 'empty'
%
% Options used (all have defaults, see bone.quality_defaults)
%   LowDensityHU, MinContrastHU, MaxLowFrac, TissueCeilingHU, MinFillRatio

if nargin < 4, opts = struct(); end
opts = bone.quality_defaults(opts);

q = struct('core_hu', NaN, 'surface_hu', NaN, 'rind_hu', NaN, ...
    'contrast', NaN, 'low_frac', NaN, 'volume_mm3', 0, ...
    'fill_ratio', 0, 'score', 0, 'flags', {{}});

if isempty(mask) || ~any(mask(:))
    q.flags = {'empty'};
    return;
end

vol = double(vol);
voxmm = mean(spacing);
q.volume_mm3 = sum(mask(:)) * prod(spacing);

% ---- Interior, surface and rind ----
D_in = bwdist(~mask) * voxmm;
interior = mask & (D_in > 1.0);
if nnz(interior) < 20
    interior = mask & (D_in > 0.5);
end
if nnz(interior) < 20
    interior = mask;
end
q.core_hu = median(vol(interior));

shell = mask & ~imerode(mask, strel('sphere', 1));
if any(shell(:))
    q.surface_hu = median(vol(shell));
end

% Rind: non-air material just beyond the mask — the tissue the segmentation
% chose to leave out. If it looks like what we kept, the boundary is not
% backed by a real density step. The shell immediately against the mask is
% skipped because partial-volume voxels there are a blend of bone and air
% and would read as tissue on any specimen.
outside = imdilate(mask, strel('sphere', 3)) & ~imdilate(mask, strel('sphere', 1));
rind = outside & (vol > -300);
if nnz(rind) >= 20
    q.rind_hu = median(vol(rind));
else
    % Nothing but air around the bone — a clean excised specimen.
    q.rind_hu = -1000;
end

q.contrast = q.core_hu - q.rind_hu;
q.low_frac = sum(vol(mask) < opts.TissueCeilingHU) / nnz(mask);

% Fraction of the blob this bone grew from that it ended up keeping
if opts.SourceVolMM3 > 0
    q.fill_ratio = q.volume_mm3 / opts.SourceVolMM3;
end

% ---- Flags ----
flags = {};
if q.core_hu < opts.LowDensityHU
    flags{end+1} = 'low_density';
end
if q.contrast < opts.MinContrastHU || q.low_frac > opts.MaxLowFrac
    flags{end+1} = 'tissue_suspect';
end
if opts.SourceVolMM3 > 0 && q.fill_ratio < opts.MinFillRatio
    flags{end+1} = 'under_segmented';
end
q.flags = flags;

% ---- Score ----
% Balanced so that growing into soft tissue (contrast down, low_frac up)
% costs more than it gains. Used only to compare candidates from the same
% scan against each other, never as an absolute quality number.
contrast_term = min(max(q.contrast, 0), 400) / 400;
density_term  = min(max(q.core_hu, 0), 600) / 600;
q.score = 0.45 * contrast_term + 0.35 * density_term - 0.40 * q.low_frac;
end
