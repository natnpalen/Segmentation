function [hu_min, hu_max] = percentile_clamp(vals_in_mask, prct)
% PERCENTILE_CLAMP  Robust min/max from values strictly inside mask.

if nargin<2 || isempty(prct), prct = [2 98]; end
vals_in_mask = vals_in_mask(isfinite(vals_in_mask));
if isempty(vals_in_mask)
    hu_min = 0; hu_max = 1; return;
end
p = prctile(vals_in_mask, prct);
hu_min = p(1); hu_max = p(2);
if ~isfinite(hu_min) || ~isfinite(hu_max) || hu_max <= hu_min
    hu_min = min(vals_in_mask); hu_max = max(vals_in_mask);
end
end
