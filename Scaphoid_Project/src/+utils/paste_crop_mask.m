function BWfull = paste_crop_mask(BWc, crop, fullSize)
% Paste cropped logical mask BWc back into full volume.
BWfull = false(fullSize);
BWfull(crop.rRange, crop.cRange, crop.sRange) = BWc;
end
