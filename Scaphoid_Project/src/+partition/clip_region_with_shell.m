function FV_out = clip_region_with_shell(FV_reg, FV_shell)
% Inputs: FV_reg, FV_shell are structs with .vertices (Nx3), .faces (Mx3), in mm
% Output: FV_out is the clipped region FV (may be empty)
%
% Speed path: cull shell triangles by AABB against region AABB before boolean.
% Padding avoids borderline exclusions; result of 'and' intersection is identical.

V1 = double(FV_reg.vertices);  F1 = double(FV_reg.faces);
V2 = double(FV_shell.vertices);F2 = double(FV_shell.faces);

% Quick cleanups help booleans a lot
[V1,F1] = meshcheckrepair(V1, F1, 'dup');
[V1,F1] = meshcheckrepair(V1, F1, 'isolated');
[V2,F2] = meshcheckrepair(V2, F2, 'dup');
[V2,F2] = meshcheckrepair(V2, F2, 'isolated');

% ---------- AABB-based ROI cull on shell (safe, output-identical) ----------
pad = 0.5; % mm slack
minR = min(V1,[],1) - pad;
maxR = max(V1,[],1) + pad;

% Shell triangle AABBs
V2F1 = V2(F2(:,1),:); V2F2 = V2(F2(:,2),:); V2F3 = V2(F2(:,3),:);
minS = min(min(V2F1, V2F2), V2F3);
maxS = max(max(V2F1, V2F2), V2F3);

% Keep only shell faces whose AABB intersects region AABB
keep = ~( maxS(:,1) < minR(1) | minS(:,1) > maxR(1) | ...
          maxS(:,2) < minR(2) | minS(:,2) > maxR(2) | ...
          maxS(:,3) < minR(3) | minS(:,3) > maxR(3) );

F2r = F2(keep,:);
if isempty(F2r)
    FV_out = struct('vertices',[], 'faces',[]);
    return;
end

% Compact shell vertices & remap faces
used = false(size(V2,1),1); used(F2r(:)) = true;
map = zeros(size(used)); map(used) = 1:nnz(used);
V2r = V2(used,:);
F2r = map(F2r);

% Optional: quick shell-only repair after compaction
[V2r,F2r] = meshcheckrepair(V2r, F2r, 'dup');
[V2r,F2r] = meshcheckrepair(V2r, F2r, 'isolated');

try
    % Primary path: exact triangle boolean (iso2mesh)
    [V,F] = surfboolean(V1, F1, 'and', V2r, F2r);   % intersection
    if isempty(V) || isempty(F)
        FV_out = struct('vertices',[], 'faces',[]);
        return;
    end

    % Final cleanup
    [V,F] = meshcheckrepair(V, F, 'dup');
    [V,F] = meshcheckrepair(V, F, 'isolated');
    [V,F] = removeisolatednode(V, F);
    FV_out = struct('vertices', V, 'faces', F);
    return;

catch ME
    % --- Fallback: inpolyhedron cull + optional snap ---
    warning('surfboolean failed (%s); using fallback.', ME.message);

    C = (V1(F1(:,1),:) + V1(F1(:,2),:) + V1(F1(:,3),:)) / 3;
    inside = inpolyhedron(F2r, V2r, C);
    Fk = F1(inside, :);
    if isempty(Fk)
        FV_out = struct('vertices',[], 'faces',[]);
        return;
    end

    try
        [~, Pclosest, sd] = point2trimesh('Faces', F2r, 'Vertices', V2r, 'QueryPoints', V1);
        near = (sd > 0) & (sd < 1.0);
        V1(near, :) = Pclosest(near, :);
    catch
        % keep culled faces as-is
    end

    [V,F] = meshcheckrepair(V1, Fk, 'dup');
    [V,F] = removeisolatednode(V, F);
    FV_out = struct('vertices', V, 'faces', F);
end
end
