function stats = label_stats(labels, HU, ds, k, mask)
% LABEL_STATS
% Compute basic per-label stats for audit & naming.
% Returns an array of structs with fields:
%   label, voxels, mm3, meanHU, medianHU
%
% Inputs
%   labels : uint* label map
%   HU     : double HU volume
%   ds     : struct with .spacing (1x3 mm)
%   k      : optional number of labels (will infer from data if omitted)
%   mask   : optional logical mask (defaults to labels>0)

if nargin < 5 || isempty(mask), mask = (labels>0); else, mask = logical(mask); end
if nargin < 4 || isempty(k), k = double(max(labels(:))); end
spacing = double(ds.spacing(:))';
voxmm3 = prod(spacing);

stats = repmat(struct('label',0,'voxels',0,'mm3',0,'meanHU',NaN,'medianHU',NaN), k, 1);

for lbl = 1:k
    msk = (labels==lbl) & mask;
    v = nnz(msk);
    stats(lbl).label  = lbl;
    stats(lbl).voxels = v;
    stats(lbl).mm3    = v * voxmm3;
    if v>0
        vals = HU(msk);
        stats(lbl).meanHU   = mean(vals, 'omitnan');
        stats(lbl).medianHU = median(vals, 'omitnan');
    end
end
end
