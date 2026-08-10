function pack_result = pack_specimens(bone_mask, cortical, cancellous, ds, stl_paths, shape_names, opts, bone_axis, bone_idx)
% PACK_SPECIMENS  Pack mechanical test specimens into bone regions.
%
%   pack_result = bone.pack_specimens(bone_mask, cortical, cancellous, ds, ...
%       stl_paths, shape_names, opts, bone_axis, bone_idx)
%
% When opts.PackWholeBone is true, packs into the full bone mask as a single
% region (ignoring cortical/cancellous boundaries). Otherwise packs cortical
% and cancellous regions separately.
%
% Uses true mesh geometry for fitting: rotates STL vertices first, then
% voxelizes the rotated mesh. Crops regions to their bounding box before
% convolution for speed, then maps results back to full-volume coordinates.

if nargin < 9, bone_idx = 0; end
tag = sprintf('Bone #%d', bone_idx);

spacing = ds.spacing;
vol = double(ds.HU);
voxel_vol = prod(spacing);
n_shapes = numel(stl_paths);

cort_vol = sum(cortical(:)) * voxel_vol;
canc_vol = sum(cancellous(:)) * voxel_vol;

% ---- Load STL meshes and build rotated templates ----
fprintf('       [%s] Building templates (%d shapes x %d orientations)...\n', ...
    tag, n_shapes, opts.PackingOrientations);
rotations = generate_bone_aligned_rotations(bone_axis, opts.PackingOrientations);
n_orient = size(rotations, 3);

templates = {};
for si = 1:n_shapes
    try
        TR = stlread(stl_paths{si});
        V_raw = double(TR.Points);
        F = double(TR.ConnectivityList);
        V_raw = V_raw - mean(V_raw, 1);
        bbox_mm = max(V_raw, [], 1) - min(V_raw, [], 1);
    catch
        continue;
    end

    % Canonical shape alignment: rotate the shape so its own longest
    % principal axis is +z BEFORE the roll/tilt/bone-frame rotation. The
    % STLs are modeled with arbitrary axes (Bend's 20mm side is X, Shear's
    % 40mm side is Z), so without this the meaning of "aligned with the
    % bone" depended on how each STL happened to be drawn.
    [Ev, ~] = eig(V_raw' * V_raw);
    s_axis = Ev(:, 3);
    if abs(s_axis(1)) < 0.9
        s_p = cross(s_axis, [1; 0; 0]);
    else
        s_p = cross(s_axis, [0; 1; 0]);
    end
    s_p = s_p / norm(s_p);
    s_p2 = cross(s_axis, s_p);
    s_p2 = s_p2 / norm(s_p2);
    S_align = [s_p, s_p2, s_axis]';   % maps shape long axis -> template +z

    shape_tpls = {};
    n_valid = 0;
    n_dup = 0;
    for oi = 1:n_orient
        R = rotations(:,:,oi) * S_align;
        V_rot = (R * V_raw')';
        [shape_mask, ~] = voxelize_mesh(V_rot, F, spacing);

        shape_vol = sum(shape_mask(:)) * voxel_vol;
        if shape_vol < 0.1, continue; end

        % Symmetry dedup: a 180-degree roll of a rectangular bar — or any
        % roll of a cylinder — voxelizes to an identical mask. Searching
        % the same voxel shape twice is pure waste, so keep one of each.
        is_dup = false;
        for di = 1:numel(shape_tpls)
            if isequal(shape_tpls{di}.sz, size(shape_mask)) && ...
                    isequal(shape_tpls{di}.mask, shape_mask)
                is_dup = true;
                break;
            end
        end
        if is_dup
            n_dup = n_dup + 1;
            continue;
        end

        tpl = struct();
        tpl.shape_idx = si;
        tpl.shape_name = shape_names{si};
        tpl.orientation = oi;
        tpl.rotation = R;
        tpl.vertices_mm = V_rot;
        tpl.faces = F;
        tpl.mask = shape_mask;
        tpl.volume_mm3 = shape_vol;
        tpl.sz = size(shape_mask);
        shape_tpls{end+1} = tpl; %#ok<AGROW>
        n_valid = n_valid + 1;
    end
    templates = [templates, shape_tpls];
    fprintf('       [%s] %s: %d unique orientations (%d tried, %d symmetric duplicates, %.0f mm3 each)\n', ...
        tag, shape_names{si}, n_valid, n_orient, n_dup, ...
        ternary(n_valid > 0, shape_tpls{end}.volume_mm3, 0));
end

if isempty(templates)
    fprintf('       [%s] No valid templates — skipping\n', tag);
    pack_result = empty_result();
    return;
end

% ---- Pack each region ----
whole_bone_mode = isfield(opts, 'PackWholeBone') && opts.PackWholeBone;

if whole_bone_mode
    region_vol = sum(bone_mask(:)) * voxel_vol;
    fprintf('       [%s] Packing whole bone (%.0f mm3)...\n', tag, region_vol);
    whole_placements = pack_region(bone_mask, templates, vol, spacing, shape_names, n_shapes, tag, 'whole');
    fprintf('       [%s] Placed %d specimens (whole bone)\n', tag, numel(whole_placements));

    pack_result = struct();
    pack_result.whole_bone = true;
    pack_result.whole_placements = whole_placements;
    pack_result.cortical_placements = empty_placements();
    pack_result.cancellous_placements = empty_placements();
    pack_result.n_whole = numel(whole_placements);
    pack_result.n_cortical = 0;
    pack_result.n_cancellous = 0;
    pack_result.n_total = numel(whole_placements);

    pack_result.summary = struct();
    for si = 1:n_shapes
        n_w = sum(arrayfun(@(p) p.shape_idx == si, whole_placements));
        pack_result.summary.(shape_names{si}) = struct('cortical', 0, 'cancellous', 0, 'whole', n_w);
    end
else
    fprintf('       [%s] Packing cortical (%.0f mm3)...\n', tag, cort_vol);
    cort_placements = pack_region(cortical, templates, vol, spacing, shape_names, n_shapes, tag, 'cortical');
    fprintf('       [%s] Packing cancellous (%.0f mm3)...\n', tag, canc_vol);
    canc_placements = pack_region(cancellous, templates, vol, spacing, shape_names, n_shapes, tag, 'cancellous');
    fprintf('       [%s] Placed %d cortical + %d cancellous\n', tag, ...
        numel(cort_placements), numel(canc_placements));

    pack_result = struct();
    pack_result.whole_bone = false;
    pack_result.whole_placements = empty_placements();
    pack_result.cortical_placements = cort_placements;
    pack_result.cancellous_placements = canc_placements;
    pack_result.n_whole = 0;
    pack_result.n_cortical = numel(cort_placements);
    pack_result.n_cancellous = numel(canc_placements);
    pack_result.n_total = numel(cort_placements) + numel(canc_placements);

    pack_result.summary = struct();
    for si = 1:n_shapes
        n_cort = sum(arrayfun(@(p) p.shape_idx == si, cort_placements));
        n_canc = sum(arrayfun(@(p) p.shape_idx == si, canc_placements));
        pack_result.summary.(shape_names{si}) = struct('cortical', n_cort, 'cancellous', n_canc, 'whole', 0);
    end
end
end


% =========================================================================
%  VOXELIZE ROTATED MESH (inline, no file I/O)
% =========================================================================
function [vox_mask, grid_sz] = voxelize_mesh(V, F, spacing)

    V_vox = zeros(size(V));
    V_vox(:,1) = V(:,1) / spacing(1);
    V_vox(:,2) = V(:,2) / spacing(2);
    V_vox(:,3) = V(:,3) / spacing(3);

    V_vox = V_vox - min(V_vox, [], 1) + 2;
    grid_sz = ceil(max(V_vox, [], 1)) + 2;

    vox_mask = false(grid_sz);
    z_min = floor(min(V_vox(:,3)));
    z_max = ceil(max(V_vox(:,3)));

    for z = max(1, z_min):min(grid_sz(3), z_max)
        z_vals = reshape(V_vox(F, 3), size(F));
        f_min = min(z_vals, [], 2);
        f_max = max(z_vals, [], 2);
        active = find(f_min <= z & f_max >= z);

        if isempty(active), continue; end

        segments = [];
        for fi = 1:numel(active)
            tri = F(active(fi), :);
            v1 = V_vox(tri(1), :);
            v2 = V_vox(tri(2), :);
            v3 = V_vox(tri(3), :);
            pts = intersect_triangle_z(v1, v2, v3, z);
            if size(pts, 1) >= 2
                segments = [segments; pts(1,:) pts(2,:)]; %#ok<AGROW>
            end
        end

        if isempty(segments), continue; end

        all_y = [segments(:,1); segments(:,3)];
        y_min_s = max(1, floor(min(all_y)));
        y_max_s = min(grid_sz(1), ceil(max(all_y)));

        for y = y_min_s:y_max_s
            x_hits = [];
            for si = 1:size(segments, 1)
                p1 = segments(si, 1:2);
                p2 = segments(si, 3:4);
                if (p1(1) <= y && p2(1) > y) || (p2(1) <= y && p1(1) > y)
                    t = (y - p1(1)) / (p2(1) - p1(1));
                    x_hit = p1(2) + t * (p2(2) - p1(2));
                    x_hits(end+1) = x_hit; %#ok<AGROW>
                end
            end
            x_hits = sort(x_hits);
            for pi = 1:2:numel(x_hits)-1
                x1 = max(1, round(x_hits(pi)));
                x2 = min(grid_sz(2), round(x_hits(pi+1)));
                if x1 <= x2
                    vox_mask(y, x1:x2, z) = true;
                end
            end
        end
    end

    if ~any(vox_mask(:))
        vox_mask = surface_fill(V_vox, F, grid_sz);
    end
end


% =========================================================================
%  REGION PACKING (with bounding-box crop for speed)
% =========================================================================
function placements = pack_region(region, templates, vol, spacing, shape_names, n_shapes, tag, region_label)

    voxel_vol = prod(spacing);
    region_vol = sum(region(:)) * voxel_vol;
    placements = struct('shape_name', {}, 'shape_idx', {}, 'orientation', {}, ...
        'position_vox', {}, 'volume_mm3', {}, 'mean_hu', {}, ...
        'vertices_mm', {}, 'faces', {});

    if region_vol < 1.0
        fprintf('       [%s] %s region too small (%.0f mm3) — skipping\n', tag, region_label, region_vol);
        return;
    end

    available = region;

    % Phase 1: one of each type (priority)
    for si = 1:n_shapes
        type_idx = find(cellfun(@(t) t.shape_idx == si, templates));
        if isempty(type_idx), continue; end

        [p, available, best_fit, skip_msg] = try_place_best(available, templates(type_idx), vol, spacing);
        if ~isempty(p)
            placements(end+1) = p; %#ok<AGROW>
            fprintf('       [%s] Placed %s in %s (%.0f mm3, HU %.0f)\n', ...
                tag, p.shape_name, region_label, p.volume_mm3, p.mean_hu);
        elseif ~isempty(skip_msg)
            fprintf('       [%s] %s does not fit in %s (%s)\n', ...
                tag, shape_names{si}, region_label, skip_msg);
        else
            fprintf('       [%s] %s does not fit in %s (best overlap %.0f%%)\n', ...
                tag, shape_names{si}, region_label, best_fit * 100);
        end
    end

    % Phase 2: greedily pack more specimens.
    % The available region only ever shrinks, so a template that has no
    % valid position now can never gain one — deactivate it permanently
    % instead of re-evaluating it every round.
    max_additional = 50;
    active = true(1, numel(templates));
    for attempt = 1:max_additional
        if ~any(available(:)) || ~any(active), break; end

        idx_active = find(active);
        [p, available, ~, ~, had_fit] = try_place_best(available, templates(idx_active), vol, spacing);
        active(idx_active(~had_fit)) = false;
        if isempty(p), break; end

        placements(end+1) = p; %#ok<AGROW>
        fprintf('       [%s] Placed %s in %s (%.0f mm3, HU %.0f) [greedy #%d]\n', ...
            tag, p.shape_name, region_label, p.volume_mm3, p.mean_hu, attempt);
    end
end


% =========================================================================
%  PLACEMENT SEARCH (convolution on cropped ROI for speed)
% =========================================================================
function [placement, available, best_fit_frac, skip_msg, had_fit] = try_place_best(available, templates, vol, spacing)

    placement = [];
    best_score = -Inf;
    best_pos = [];
    best_tpl = [];
    best_fit_frac = 0;
    skip_msg = '';
    had_fit = false(1, numel(templates));

    full_sz = size(available);

    % Crop available region to its bounding box (with padding for templates).
    % A valid placement needs >=95% of the specimen inside the available
    % region, so a template can only overhang the bounding box by a small
    % fraction of its own size — a quarter-size pad is sufficient.
    max_tpl_sz = [0 0 0];
    for ti = 1:numel(templates)
        max_tpl_sz = max(max_tpl_sz, templates{ti}.sz);
    end

    [rr, cc, ss] = ind2sub(full_sz, find(available));
    if isempty(rr), return; end
    pad = max(4, ceil(max_tpl_sz * 0.25));
    roi_min = max([1 1 1], [min(rr) min(cc) min(ss)] - pad);
    roi_max = min(full_sz, [max(rr) max(cc) max(ss)] + pad);

    avail_crop = available(roi_min(1):roi_max(1), roi_min(2):roi_max(2), roi_min(3):roi_max(3));
    vol_crop = vol(roi_min(1):roi_max(1), roi_min(2):roi_max(2), roi_min(3):roi_max(3));
    crop_sz = size(avail_crop);

    D_crop = bwdist(~avail_crop) .* mean(spacing);

    avail_vol_mm3 = sum(avail_crop(:)) * prod(spacing);

    % FFT-based placement search. The overlap count at every position is a
    % correlation of the availability map with the template mask; doing it
    % in the frequency domain costs O(N log N) regardless of template size,
    % where spatial convn costs O(N * template_vox) — minutes per attempt
    % for large templates like Shear. The availability FFT is shared across
    % all templates; the depth-map FFT is computed lazily on first use.
    fftsz = crop_sz + max_tpl_sz - 1;
    FA = fftn(single(avail_crop), fftsz);
    FD = [];

    n_eval = 0;
    n_skip_vol = 0;
    n_skip_bbox = 0;
    min_skip_vol = Inf;

    for ti = 1:numel(templates)
        tpl = templates{ti};
        tsz = tpl.sz;

        if any(tsz > crop_sz)
            n_skip_bbox = n_skip_bbox + 1;
            continue;
        end

        % Quick volume check: skip if template is larger than remaining region
        if tpl.volume_mm3 > avail_vol_mm3 * 1.1
            n_skip_vol = n_skip_vol + 1;
            min_skip_vol = min(min_skip_vol, tpl.volume_mm3);
            continue;
        end
        n_eval = n_eval + 1;

        n_template_vox = sum(tpl.mask(:));
        FK = fftn(flip_3d(single(tpl.mask)), fftsz);
        overlap_full = real(ifftn(FA .* FK));
        overlap_count = overlap_full(tsz(1):crop_sz(1), tsz(2):crop_sz(2), tsz(3):crop_sz(3));
        % FFT round-trip noise is well under one voxel count — snap to
        % integers so the fit fraction is exact.
        overlap_count = max(0, round(overlap_count));
        fit_frac = overlap_count / max(1, n_template_vox);

        % Track best fit fraction across all templates (for diagnostics)
        max_fit_this = max(fit_frac(:));
        if max_fit_this > best_fit_frac
            best_fit_frac = max_fit_this;
        end

        good_fit = fit_frac >= 0.95;
        if ~any(good_fit(:)), continue; end
        had_fit(ti) = true;

        if isempty(FD)
            FD = fftn(single(D_crop), fftsz);
        end
        depth_full = real(ifftn(FD .* FK));
        depth_sum = depth_full(tsz(1):crop_sz(1), tsz(2):crop_sz(2), tsz(3):crop_sz(3));
        avg_depth = depth_sum / max(1, n_template_vox);

        score_map = fit_frac + avg_depth;
        score_map(~good_fit) = -Inf;

        [local_best, linear_idx] = max(score_map(:));
        if local_best > best_score
            [pr, pc, ps] = ind2sub(size(score_map), linear_idx);
            best_score = local_best;
            best_pos = [pr + roi_min(1) - 1, pc + roi_min(2) - 1, ps + roi_min(3) - 1];
            best_tpl = tpl;
        end
    end

    if isempty(best_pos)
        % Explain why nothing was even evaluated, so "0% overlap" is never
        % printed for a specimen that was skipped before evaluation.
        if n_eval == 0 && n_skip_vol > 0
            skip_msg = sprintf('specimen volume %.0f mm3 exceeds available region %.0f mm3', ...
                min_skip_vol, avail_vol_mm3);
        elseif n_eval == 0 && n_skip_bbox > 0
            skip_msg = 'specimen extent exceeds region bounding box';
        end
        return;
    end

    r1 = best_pos(1); c1 = best_pos(2); s1 = best_pos(3);
    tsz = best_tpl.sz;
    r2 = r1 + tsz(1) - 1;
    c2 = c1 + tsz(2) - 1;
    s2 = s1 + tsz(3) - 1;

    placed_mask = false(full_sz);
    placed_mask(r1:r2, c1:c2, s1:s2) = best_tpl.mask;
    placed_mask = placed_mask & available;

    if sum(placed_mask(:)) < 0.90 * sum(best_tpl.mask(:))
        return;
    end

    % Reconstruct placed mesh vertices in volume mm coordinates.
    % Voxelization does: V_vox = V_mm./spacing - min(V_mm./spacing) + 2
    % So template voxel = V_mm./spacing - min(V_mm./spacing) + 2
    % Template voxel (i,j,k) placed at volume voxel (r1+i-1, c1+j-1, s1+k-1)
    % Volume voxel of vertex = V_mm./spacing - min(V_mm./spacing) + 2 + [r1-1, c1-1, s1-1]
    % Volume mm = volume_voxel .* spacing
    V_mm = best_tpl.vertices_mm;
    V_vox_local = [V_mm(:,1)/spacing(1), V_mm(:,2)/spacing(2), V_mm(:,3)/spacing(3)];
    min_vox = min(V_vox_local, [], 1);
    V_vox_template = V_vox_local - min_vox + 2;
    V_vol_vox = V_vox_template + [r1-1, c1-1, s1-1];
    V_placed = [V_vol_vox(:,1)*spacing(1), V_vol_vox(:,2)*spacing(2), V_vol_vox(:,3)*spacing(3)];

    placement = struct();
    placement.shape_name = best_tpl.shape_name;
    placement.shape_idx = best_tpl.shape_idx;
    placement.orientation = best_tpl.orientation;
    placement.position_vox = best_pos;
    placement.volume_mm3 = best_tpl.volume_mm3;
    placement.mean_hu = mean(vol(placed_mask));
    placement.vertices_mm = V_placed;
    placement.faces = best_tpl.faces;

    available(placed_mask) = false;
end


% =========================================================================
%  HELPER FUNCTIONS
% =========================================================================
function flipped = flip_3d(A)
    flipped = flip(flip(flip(A, 1), 2), 3);
end


function pts = intersect_triangle_z(v1, v2, v3, z)
    edges = {v1, v2; v2, v3; v3, v1};
    pts = zeros(0, 2);
    for e = 1:3
        p1 = edges{e, 1};
        p2 = edges{e, 2};
        if (p1(3) <= z && p2(3) >= z) || (p2(3) <= z && p1(3) >= z)
            dz = p2(3) - p1(3);
            if abs(dz) < 1e-10
                pts = [pts; p1(1:2); p2(1:2)]; %#ok<AGROW>
            else
                t = (z - p1(3)) / dz;
                t = max(0, min(1, t));
                pt = p1 + t * (p2 - p1);
                pts = [pts; pt(1:2)]; %#ok<AGROW>
            end
        end
    end
end


function mask = surface_fill(V, F, grid_sz)
    mask = false(grid_sz);
    for fi = 1:size(F, 1)
        v1 = V(F(fi,1), :);
        v2 = V(F(fi,2), :);
        v3 = V(F(fi,3), :);
        for u = 0:0.1:1
            for w = 0:0.1:(1-u)
                pt = v1*(1-u-w) + v2*u + v3*w;
                idx = round(pt);
                if all(idx >= 1) && idx(1) <= grid_sz(1) && ...
                   idx(2) <= grid_sz(2) && idx(3) <= grid_sz(3)
                    mask(idx(1), idx(2), idx(3)) = true;
                end
            end
        end
    end
    mask = imfill(mask, 'holes');
end


% =========================================================================
%  BONE-ALIGNED ROTATIONS
% =========================================================================
function R = generate_bone_aligned_rotations(bone_axis, n_orient)

    bone_axis = bone_axis(:) / norm(bone_axis);

    if abs(bone_axis(1)) < 0.9
        perp = cross(bone_axis, [1;0;0]);
    else
        perp = cross(bone_axis, [0;1;0]);
    end
    perp = perp / norm(perp);
    perp2 = cross(bone_axis, perp);
    perp2 = perp2 / norm(perp2);

    % Columns [perp, perp2, bone_axis]: maps template +z onto the bone's
    % long axis, so phi=0 orientations lie ALONG the bone and theta rolls
    % about it. The previous transposed form ([...]') mapped the bone axis
    % onto +z instead — template orientations relative to the bone were
    % scrambled, and with 6 orientations only one accidentally aligned an
    % elongated specimen with the shaft (the likely reason Bend and Shear
    % never fit whole metacarpals).
    B_bone = [perp, perp2, bone_axis];

    % Ordered (roll about bone axis, tilt off axis) pairs. The first 8
    % reproduce the original coarse set (90-degree steps), so low
    % PackingOrientations values behave exactly as before. The pairs after
    % that add finer rolls and small tilts: for an elongated specimen in a
    % curved bone, the difference between fitting and not fitting is often
    % a 15-45 degree roll or tilt that 90-degree steps can never reach.
    pairs = [ ...
          0  0;   0 90;  90  0;  90 90; 180  0; 180 90; 270  0; 270 90; ... % original coarse set
         30  0;  45  0;  60  0; 120  0; 135  0; 150  0; ...                 % finer rolls
        210  0; 225  0; 240  0; 300  0; 315  0; 330  0; ...
          0 15;  90 15; 180 15; 270 15; ...                                 % small tilts
          0 -15; 90 -15; 180 -15; 270 -15; ...
          0 75;  90 75; 180 75; 270 75; ...                                 % near-perpendicular tilts
          0 105; 90 105];

    n = max(1, min(n_orient, size(pairs, 1)));
    R = zeros(3, 3, n);

    for idx = 1:n
        theta = deg2rad(pairs(idx, 1));
        phi = deg2rad(pairs(idx, 2));

        Rz = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];

        R(:,:,idx) = B_bone * Rz * Rx;
    end
end


% =========================================================================
%  UTILITIES
% =========================================================================
function r = empty_result()
    r = struct();
    r.whole_bone = false;
    r.whole_placements = empty_placements();
    r.cortical_placements = empty_placements();
    r.cancellous_placements = empty_placements();
    r.n_whole = 0;
    r.n_cortical = 0;
    r.n_cancellous = 0;
    r.n_total = 0;
    r.summary = struct();
end

function p = empty_placements()
    p = struct('shape_name', {}, 'shape_idx', {}, ...
        'orientation', {}, 'position_vox', {}, 'volume_mm3', {}, ...
        'mean_hu', {}, 'vertices_mm', {}, 'faces', {});
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
