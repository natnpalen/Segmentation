function D = distanceToLabelVox(labels, targetLabel)
% Euclidean voxel distance to the target label (0 means on/inside target).
% Returns double array in *voxel* units (we only need equality of distances).
LBL = uint16(labels);
inside = (LBL == uint16(targetLabel));
if ~any(inside(:))
    D = inf(size(LBL), 'double');
    return;
end
% distance *to* target label: inside has 0, others get distance to nearest target voxel
D = bwdist(inside, 'euclidean');   % voxel-units (exact EDT on grid)
end
