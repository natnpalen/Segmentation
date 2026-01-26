function Vc = apply_crop(V, crop)
% Crop any 3-D volume to crop ranges.
Vc = V(crop.rRange, crop.cRange, crop.sRange);
end