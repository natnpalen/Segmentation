function grad = normalize_hu_to_density(HU, mask, hu_min, hu_max)
% NORMALIZE_HU_TO_DENSITY  Map HU → [0,1] inside mask; NaN outside.
rng = max(hu_max - hu_min, eps);
grad = nan(size(HU), 'double');
inb  = mask & isfinite(HU);
grad(inb) = (double(HU(inb)) - hu_min) / rng;
grad(inb) = min(max(grad(inb), 0), 1);
end
