function sdf = signedDistanceMM(BW, spacing)
% signedDistanceMM: positive inside (foreground), negative outside, in mm.
% Uses exact Euclidean distance transforms in voxels, scaled to mm.
BW = logical(BW);
d_in  = bwdist(~BW);   % distance from inside voxels to boundary (vox)
d_out = bwdist(BW);    % distance from outside to boundary (vox)
voxel_mm = mean(spacing);  % isotropic approx is fine for surface location
sdf = (d_in - d_out) * voxel_mm;
end
