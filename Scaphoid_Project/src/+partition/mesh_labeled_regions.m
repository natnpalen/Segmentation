function meshes = mesh_labeled_regions(label_volume, ds_hires, k, HU_hires, ~)
% MESH_LABELED_REGIONS — Per-label watertight meshes using
% "closest-other + mask-complement" distance formulation (no outer SDF).
%
% INPUTS
%   label_volume : uint* 3D label array (0 = air)
%   ds_hires.spacing : 1x3 mm
%   k : number of labels (1..k)
%   HU_hires : unused (kept for signature compatibility)
% OUTPUT
%   meshes : kx1 cell of structs {vertices (mm), faces}

spacing = ds_hires.spacing(:)'; % [sx sy sz] mm
meshes = cell(k, 1);

L = uint16(label_volume);
PAIR_BAND_VOX = 3; % narrow band around interfaces (voxels)

% Legacy outer mask for speed (anything >0 considered inside)
outer_mask = label_volume > 0;

% ---------- Distances (voxel units) ----------
Dlabel = cell(k+1, 1); % 1..k are labels, k+1 is 0 (air)
parfor lab = 1:k
    Dlabel{lab} = meshing.distanceToLabelVox(L, lab);
end
Dlabel{k+1} = meshing.distanceToLabelVox(L, 0);

% Distance to mask complement (outside of scaphoid)
outside_mask = ~outer_mask;
D_outside = meshing.distanceToLabelVox(uint16(outside_mask), 1);
D_outside_d = double(D_outside); % small pre-cast once

% ================== Per-label extraction ==================
parfor lab = 1:k   % <<< parallelized (output order preserved via indexing)
    Di = double(Dlabel{lab});

    % Min distance to competitors: other labels + air + outside mask
    Dmin_other = D_outside_d; % start with mask competitor
    for m = [1:lab-1, lab+1:k+1]
        Dm = double(Dlabel{m});
        % elementwise min to build "closest other"
        Dmin_other = min(Dmin_other, Dm);
    end

    % Implicit field: positive outside ℓ, negative inside ℓ
    S = Dmin_other - Di;

    % Narrow band mask for speed/stability
    Near = (min(Di, Dmin_other) <= PAIR_BAND_VOX);

    % Strict sign-change guard
    Sloc = S(Near); Sloc = Sloc(isfinite(Sloc));
    if isempty(Sloc) || ~(any(Sloc < 0) && any(Sloc > 0))
        meshes{lab} = struct('vertices',[],'faces',[]);
        continue;
    end

    % Extract surface in mm using your meshing utility
    fv = meshing.isosurface_mm(S, 0.0, spacing);
    if isempty(fv) || isempty(fv.vertices) || isempty(fv.faces)
        meshes{lab} = struct('vertices',[],'faces',[]);
        continue;
    end

    V = fv.vertices;
    F = fv.faces;

    % ---------- Mesh hardening ----------
    % Drop NaN/Inf vertices & faces that reference them
    badV = any(~isfinite(V), 2);
    if any(badV)
        keepV = ~badV;
        map = zeros(size(keepV));
        map(keepV) = 1:nnz(keepV);
        V = V(keepV, :);
        F = map_safe_faces(F, map);
    end

    % Global vertex merge with tolerance
    [V_merged, ~, remap] = uniquetol(V, 1e-5, 'ByRows', true);
    F = remap(F);

    % Drop exactly degenerate faces (repeated indices)
    F = F( F(:,1)~=F(:,2) & F(:,1)~=F(:,3) & F(:,2)~=F(:,3), : );

    % Orientation-agnostic duplicate face removal
    F = unique_faces_unoriented(F);

    % Cull tiny sliver faces by geometric area (2*area threshold)
    eps_area2 = 1e-9; % units of (2*area) in mm^2
    if ~isempty(F)
        keep = face_area2_mask(V_merged, F, eps_area2);
        F = F(keep, :);
    end

    % Ensure consistent outward winding (positive signed volume)
    if ~isempty(F) && signed_volume(V_merged, F) < 0
        F = F(:, [1 3 2]);
    end

    meshes{lab} = struct('vertices', V_merged, 'faces', F);
end

end % function mesh_labeled_regions

% =========================
% ===== Helper Utils ======
% =========================

function F2 = map_safe_faces(F, map)
% Remap faces through 'map' (zeros are invalid). Drops any face that
% references an invalid vertex after remapping.
Fm = map(F);
bad = any(Fm==0, 2);
F2 = Fm(~bad, :);
end

function F_unique = unique_faces_unoriented(F)
% Remove duplicate triangles regardless of vertex order.
if isempty(F), F_unique = F; return; end
F_sorted = sort(F, 2);
[~, ia] = unique(F_sorted, 'rows', 'stable');
F_unique = F(ia, :);
end

function keep = face_area2_mask(V, F, eps_area2)
% Keep mask for faces with 2*area > eps_area2 (area2 is norm of cross product)
if isempty(F)
    keep = false(0,1);
    return;
end
A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
area2 = vecnorm(cross(B-A, C-A, 2), 2, 2); % equals 2*area
keep = area2 > eps_area2 & isfinite(area2);
end

function sVol = signed_volume(V, F)
% Positive if triangles are outward-facing for a closed mesh.
if isempty(F)
    sVol = 0.0;
    return;
end
A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
sVol = sum(dot(A, cross(B, C, 2), 2)) / 6.0;
end
