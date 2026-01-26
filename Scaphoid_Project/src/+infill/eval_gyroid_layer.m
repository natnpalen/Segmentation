function [C, xv, yv, Phi] = eval_gyroid_layer(Z_mm, rhoGrid, origin, spacing, w, rho_eps, xy_step)
% Evaluate φ(x,y,Z) on a grid and extract φ=0 contours (marching squares).
% Returns MATLAB contour matrix C plus the sampled grid and the φ field.

% X/Y ranges from the density grid bounding box
sz = size(rhoGrid);
xmin = origin(1); xmax = origin(1) + spacing(1)*(sz(1)-1);
ymin = origin(2); ymax = origin(2) + spacing(2)*(sz(2)-1);

xv = xmin:xy_step:xmax;
yv = ymin:xy_step:ymax;

% Start a pool if not present
if isempty(gcp('nocreate'))
    try, parpool; catch, end
end

Phi = zeros(numel(yv), numel(xv));
% parfor over rows (sliced assignment => parfor-friendly)
parfor j = 1:numel(yv)
    row = zeros(1, numel(xv));
    y = yv(j);
    for i = 1:numel(xv)
        x = xv(i);
        p = [x, y, Z_mm];
        row(i) = infill.gyroid_phi_trilinear(p, rhoGrid, origin, spacing, w, rho_eps);
    end
    Phi(j,:) = row;
end

% φ = 0 contours (compute only; do not draw)
C = contourc(xv, yv, Phi, [0 0]);
end
