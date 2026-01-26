function win = compute_spacing_window(neigh_mm, spacing)
% COMPUTE_SPACING_WINDOW  Odd voxel window corresponding to neigh_mm in each axis.
w = max(1, round( double(neigh_mm) ./ double(spacing(:)') )); % [wx wy wz]
win = 2*floor(w/2) + 1;  % force odd
end
