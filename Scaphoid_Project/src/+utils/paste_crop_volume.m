function Vfull = paste_crop_volume(Vc, crop, fullSize)
Vfull = zeros(fullSize, 'like', Vc);
Vfull(crop.rRange, crop.cRange, crop.sRange) = Vc;
end