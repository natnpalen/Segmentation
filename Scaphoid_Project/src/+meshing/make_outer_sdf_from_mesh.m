function phi_mm = make_outer_sdf_from_mesh(mesh_outer, gridSize, ds, mask, band_mm)
% MAKE_OUTER_SDF_FROM_MESH  Build a *smooth* signed distance field (mm)
% to the outer mesh directly from triangle geometry (no raster EDT).
%
% Inputs
%   mesh_outer : struct with fields .vertices (Nx3, mm), .faces (Mx3)
%   gridSize   : [R C S] of target grid
%   ds         : struct with field spacing (1x3, mm)
%   mask       : logical 3D mask on the same grid (inside=true)
%   band_mm    : build φ reliably within +/- band_mm of the surface
%
% Output
%   phi_mm     : 3D single, signed distance in mm (φ<0 inside mask)
%
% Notes
% - We compute *Euclidean* distance to triangles (analytic), not EDT.
% - For speed, we only evaluate voxels whose centers are within band_mm
%   of the mesh AABB expanded by band_mm, then do a centroid k-NN prefilter.

spacing = ds.spacing(:)'; % [sx sy sz] mm
R = gridSize(1); C = gridSize(2); S = gridSize(3);

V = double(mesh_outer.vertices);  % mm
F = double(mesh_outer.faces);

% ---- Build triangle data ----
A = V(F(:,1),:); B = V(F(:,2),:); Cc = V(F(:,3),:);
TriCentroids = (A + B + Cc) / 3;
% Axis-aligned bounding box (mm)
mins = min(V,[],1) - band_mm; 
maxs = max(V,[],1) + band_mm;

% ---- Prepare voxel centers (mm) only where needed ----
% Compute index ranges overlapping the expanded AABB
ir = max(1, floor(mins(1)/spacing(1))+1) : min(R, ceil(maxs(1)/spacing(1))+1);
jc = max(1, floor(mins(2)/spacing(2))+1) : min(C, ceil(maxs(2)/spacing(2))+1);
ks = max(1, floor(mins(3)/spacing(3))+1) : min(S, ceil(maxs(3)/spacing(3))+1);

% Generate voxel center coordinates (mm) for this window
[I,J,K] = ndgrid(ir, jc, ks);
P = [ (I-1)*spacing(1), (J-1)*spacing(2), (K-1)*spacing(3) ];  % (#win vox) x 3

% ---- Nearest triangle prefilter (k-NN on centroids) ----
% Build KD-tree on centroids
Mdl = createns(TriCentroids, 'NSMethod','kdtree');
kNear = 24;  % candidate triangles per query (tune if needed)

% Batch to limit memory
batch = 200000; n = size(P,1);
dist_abs = inf(n,1);

for a = 1:batch:n
    b = min(n, a+batch-1);
    Pb = P(a:b,:);
    idx = knnsearch(Mdl, Pb, 'K', kNear);  % (b-a+1) x kNear

    % Compute exact point->triangle distances against the candidate set
    dmin = inf(b-a+1,1);
    for col = 1:size(idx,2)
        t = idx(:,col);
        % triangles for this column
        At = A(t,:); Bt = B(t,:); Ct = Cc(t,:);
        d = point_triangle_distance(Pb, At, Bt, Ct); % Euclidean
        dmin = min(dmin, d);
    end
    dist_abs(a:b) = dmin;
end

% Clamp to band (outside of window we fill with +/- band_mm)
dist_abs_clamped = min(dist_abs, band_mm);

% ---- Write distances into full volume & assign signs from mask ----
phi_mm = single( band_mm * ones(R,C,S, 'single') ); % default +band (outside)
% Fill window
phi_sub = reshape(dist_abs_clamped, numel(ir), numel(jc), numel(ks));
phi_mm(ir, jc, ks) = single(phi_sub);

% Sign from mask: inside => negative
inside = false(R,C,S); inside(ir, jc, ks) = mask(ir, jc, ks);
phi_mm(inside) = -phi_mm(inside);

end

% ======= Geometry helper: Euclidean point-triangle distance (vectorized) =======
function d = point_triangle_distance(P, A, B, C)
% Compute Euclidean distance from each point in P to its corresponding triangle ABC.
% P, A, B, C: (N x 3)
% Returns d: (N x 1)
%
% Based on "Real-Time Collision Detection" (Christer Ericson), vectorized.

AB = B - A; AC = C - A; AP = P - A;
d1 = dot(AB, AP, 2);
d2 = dot(AC, AP, 2);
% Check vertex region outside A
isA = (d1 <= 0) & (d2 <= 0);
distA = sqrt(sum((P - A).^2,2));

BP = P - B;
d3 = dot(AB, BP, 2);
d4 = dot(AC, BP, 2);
% Check vertex region outside B
isB = (d3 >= 0) & (d4 - d3 <= 0);
distB = sqrt(sum((P - B).^2,2));

CP = P - C;
d5 = dot(AB, CP, 2);
d6 = dot(AC, CP, 2);
% Check vertex region outside C
isC = (d6 >= 0) & (d5 - d6 >= 0);
distC = sqrt(sum((P - C).^2,2));

% Edge regions
vc = d1.*d4 - d3.*d2;
isAB = (vc <= 0) & (d1 >= 0) & (d3 <= 0);
v  = d1 ./ (d1 - d3 + eps);
projAB = A + v .* AB;
distAB = sqrt(sum((P - projAB).^2,2));

vb = d5.*d2 - d1.*d6;
isAC = (vb <= 0) & (d2 >= 0) & (d6 <= 0);
w  = d2 ./ (d2 - d6 + eps);
projAC = A + w .* AC;
distAC = sqrt(sum((P - projAC).^2,2));

va = d3.*d6 - d5.*d4;
isBC = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0);
w2 = (d4 - d3) ./ ((d4 - d3) + (d5 - d6) + eps);
projBC = B + w2 .* (C - B);
distBC = sqrt(sum((P - projBC).^2,2));

% Face region
isFace = ~(isA | isB | isC | isAB | isAC | isBC);
% Project onto plane and compute perp distance
n = cross(AB, AC, 2);
n_norm = sqrt(sum(n.^2,2)) + eps;
distFace = abs(dot(P - A, n, 2)) ./ n_norm;

% combine regions
d = zeros(size(distA));
d(isA) = distA(isA);
d(isB) = distB(isB);
d(isC) = distC(isC);
d(isAB) = distAB(isAB);
d(isAC) = distAC(isAC);
d(isBC) = distBC(isBC);
d(isFace) = distFace(isFace);
end
