function result = separate_bones(ds, opts)
% SEPARATE_BONES  Isolate individual bones from a multi-bone excised-in-air CT scan.
%
%   result = bone.separate_bones(ds, opts)
%
% Markers are large structured radiographic assemblies (metal letter + tabs
% + housing block, up to ~2000 mm^3).  Strategy:
%   1. Detect marker assemblies: metal (HU>1200) + surrounding fixture
%   2. Build exclusion zone: metal + physical 2mm buffer to capture housing
%   3. Find bone components in non-air minus exclusion zone (no closing)
%   4. Classify: bones vs residual fixture fragments
%   5. Per-bone refinement with hole filling

vol = double(ds.HU);
spacing = ds.spacing;
voxel_vol = prod(spacing);

% ---- Stage 1: Marker assembly detection ----
fprintf('  [Separate] Stage 1: Marker assembly detection...\n');
metal = vol > opts.TagHUMin;
[marker_mask, artifact_w] = marker_and_artifact_maps(vol, opts.MarkerRangeHU, ...
    opts.ArtifactSigmaMM, spacing);

% Marker assembly = metal + collar (HU 200-700 near metal)
% + physical 2mm buffer to capture non-metallic fixture/housing
assembly_excl = physical_dilate(metal, 2.0, spacing);

n_metal = sum(metal(:));
n_excl = sum(assembly_excl(:));
fprintf('    Metal voxels (HU>%d): %d (%.1f mm^3)\n', opts.TagHUMin, n_metal, n_metal*voxel_vol);
fprintf('    Assembly exclusion (metal + 2mm): %d voxels (%.1f mm^3)\n', n_excl, n_excl*voxel_vol);

% Identify individual marker assemblies (connected components of metal)
CC_metal = bwconncomp(metal, 26);
min_tag_vox = max(5, round(2.0 / voxel_vol));
real_tags = {};
for i = 1:CC_metal.NumObjects
    if numel(CC_metal.PixelIdxList{i}) >= min_tag_vox
        [rr, cc, ss] = ind2sub(size(vol), CC_metal.PixelIdxList{i});
        tag = struct();
        tag.label = numel(real_tags) + 1;
        tag.centroid_mm = [mean(rr) mean(cc) mean(ss)] .* spacing;
        tag.volume_mm3 = numel(CC_metal.PixelIdxList{i}) * voxel_vol;
        real_tags{end+1} = tag; %#ok<AGROW>
    end
end
fprintf('    Marker assemblies: %d\n', numel(real_tags));
for t = 1:numel(real_tags)
    fprintf('      Marker %d: %.1f mm^3 metal at [%.1f %.1f %.1f] mm\n', ...
        t, real_tags{t}.volume_mm3, real_tags{t}.centroid_mm);
end

% ---- Stage 2: Find bone candidates ----
fprintf('  [Separate] Stage 2: Finding bone candidates...\n');

% Non-air WITHOUT closing (bones are solid objects separated by air,
% closing merges nearby markers with bones)
non_air = vol > -500;
fprintf('    Non-air voxels: %d (%.0f mm^3)\n', sum(non_air(:)), sum(non_air(:))*voxel_vol);

% Subtract marker exclusion zone
bone_candidates = non_air & ~assembly_excl;
fprintf('    After marker exclusion: %d voxels (%.0f mm^3)\n', ...
    sum(bone_candidates(:)), sum(bone_candidates(:))*voxel_vol);

% Connected components
CC = bwconncomp(bone_candidates, 26);
comp_vols = cellfun(@numel, CC.PixelIdxList) * voxel_vol;
[~, vol_order] = sort(comp_vols, 'descend');
fprintf('    Connected components: %d (largest: %.0f mm^3)\n', ...
    CC.NumObjects, max(comp_vols));

% ---- Stage 3: Classify and refine ----
fprintf('  [Separate] Stage 3: Classifying and refining...\n');
bones = {};
small_count = 0;
fixture_count = 0;

% Build metal distance map for fixture detection
d_metal_mm = bwdist(metal) .* mean(spacing);

for ii = 1:CC.NumObjects
    i = vol_order(ii);
    comp = false(size(vol));
    comp(CC.PixelIdxList{i}) = true;
    comp_vol = comp_vols(i);

    if comp_vol < opts.MinBoneVolMM3
        small_count = small_count + 1;
        continue;
    end

    % HU statistics
    hu_vals = vol(comp);
    n_vox = numel(hu_vals);
    mean_hu_all = mean(hu_vals);

    % Bone-tissue HU (excluding air trapped inside)
    tissue_vals = hu_vals(hu_vals > -200);
    if ~isempty(tissue_vals)
        mean_hu_tissue = mean(tissue_vals);
        tissue_frac = numel(tissue_vals) / n_vox;
    else
        mean_hu_tissue = mean_hu_all;
        tissue_frac = 0;
    end

    % Dense bone fraction (> 200 HU)
    dense_frac = sum(hu_vals > 200) / n_vox;

    % Distance to nearest metal: check if this component is a fixture fragment
    % sitting right at the edge of the exclusion zone
    comp_dists = d_metal_mm(comp);
    median_metal_dist = median(comp_dists);
    min_metal_dist = min(comp_dists);

    fprintf('    Component %d: %.0f mm^3, mean HU %.0f (tissue %.0f)\n', ...
        i, comp_vol, mean_hu_all, mean_hu_tissue);
    fprintf('      tissue=%.0f%%, dense=%.0f%%, dist_to_metal=[%.1f, med %.1f] mm\n', ...
        tissue_frac*100, dense_frac*100, min_metal_dist, median_metal_dist);

    % Fixture fragment detection:
    % Residual fixture housing sits right at the edge of the 2mm exclusion zone.
    % It has: small volume, very close to metal, low tissue fraction
    if comp_vol < 1500 && median_metal_dist < 4.0 && tissue_frac < 0.6
        fprintf('      -> FIXTURE FRAGMENT (small, near metal, low tissue)\n');
        fixture_count = fixture_count + 1;
        continue;
    end

    % Must have some actual bone tissue
    if dense_frac < 0.02
        fprintf('      -> NOISE (%.1f%% dense)\n', dense_frac*100);
        continue;
    end

    fprintf('      -> BONE CANDIDATE\n');

    % ---- Per-bone refinement ----
    refined = comp;

    % Small closing (0.5mm) to seal surface porosity only
    refined = physical_close(refined, 0.5, spacing);

    % Per-slice 2D fill to capture enclosed marrow cavities
    for z = 1:size(refined, 3)
        sl = refined(:,:,z);
        if any(sl(:))
            refined(:,:,z) = imfill(sl, 'holes');
        end
    end
    refined = imfill(refined, 'holes');

    % Clip to stay outside marker exclusion zone
    refined = refined & ~assembly_excl;

    % Remove any HU > 1200 voxels (residual metal/artifact in bone)
    refined = refined & (vol <= opts.TagHUMin);

    % Keep only the largest connected piece
    CC_ref = bwconncomp(refined, 26);
    if CC_ref.NumObjects > 1
        ref_sizes = cellfun(@numel, CC_ref.PixelIdxList);
        [~, largest] = max(ref_sizes);
        refined = false(size(vol));
        refined(CC_ref.PixelIdxList{largest}) = true;
        fprintf('      Kept largest of %d fragments (dropped %.0f mm^3)\n', ...
            CC_ref.NumObjects, (sum(ref_sizes) - ref_sizes(largest))*voxel_vol);
    end

    refined_vol = sum(refined(:)) * voxel_vol;
    if refined_vol < opts.MinBoneVolMM3
        fprintf('      -> Too small after refinement (%.0f mm^3)\n', refined_vol);
        continue;
    end

    % HU of refined bone (tissue voxels only)
    bone_vals = vol(refined & (vol > -200));
    if ~isempty(bone_vals)
        refined_hu = mean(bone_vals);
    else
        refined_hu = mean(vol(refined));
    end

    % Centroid
    [rr, cc, ss] = ind2sub(size(refined), find(refined));
    centroid_mm = [mean(rr), mean(cc), mean(ss)] .* spacing;

    % Bounding box
    bbox = [min(rr) min(cc) min(ss) max(rr) max(cc) max(ss)];

    bone_info = struct();
    bone_info.mask = refined;
    bone_info.label = i;
    bone_info.centroid_mm = centroid_mm;
    bone_info.volume_mm3 = refined_vol;
    bone_info.mean_hu = refined_hu;
    bone_info.dense_fraction = dense_frac;
    bone_info.bbox = bbox;
    bone_info.tag_id = [];
    bone_info.tag_dist = Inf;

    bones{end+1} = bone_info; %#ok<AGROW>

    fprintf('      Final: %.0f -> %.0f mm^3, mean HU %.0f\n', ...
        comp_vol, refined_vol, refined_hu);
end

fprintf('    Filtered: %d small (< %.0f mm^3), %d fixture fragments\n', ...
    small_count, opts.MinBoneVolMM3, fixture_count);

% ---- Stage 4: Tag association ----
fprintf('  [Separate] Stage 4: Tag association...\n');
bones = associate_tags(bones, real_tags);

% Sort by volume (largest first)
vols = cellfun(@(b) b.volume_mm3, bones);
[~, order] = sort(vols, 'descend');
bones = bones(order);

fprintf('\n  Found %d bones and %d markers in scan\n', numel(bones), numel(real_tags));
for i = 1:numel(bones)
    b = bones{i};
    if ~isempty(b.tag_id)
        tag_str = sprintf('marker %d (%.1f mm away)', b.tag_id, b.tag_dist);
    else
        tag_str = 'no marker';
    end
    fprintf('    Bone %d: %.1f mm^3, mean HU %.0f, centroid [%.1f %.1f %.1f] mm, %s\n', ...
        i, b.volume_mm3, b.mean_hu, b.centroid_mm, tag_str);
end

% Build specimen mask (union of all bones)
specimen = false(size(vol));
for i = 1:numel(bones)
    specimen = specimen | bones{i}.mask;
end

result = struct();
result.bones = bones;
result.specimen = specimen;
result.marker_mask = marker_mask;
result.artifact_weight = artifact_w;
result.n_tags = numel(real_tags);
end


% =========================================================================
function dilated = physical_dilate(mask, radius_mm, spacing)
    radius_vox = ceil(radius_mm ./ spacing);
    [Y, X, Z] = ndgrid(-radius_vox(1):radius_vox(1), ...
                        -radius_vox(2):radius_vox(2), ...
                        -radius_vox(3):radius_vox(3));
    dist_mm_sq = (Y*spacing(1)).^2 + (X*spacing(2)).^2 + (Z*spacing(3)).^2;
    se = strel(dist_mm_sq <= radius_mm^2);
    dilated = imdilate(mask, se);
end


function closed = physical_close(mask, radius_mm, spacing)
    radius_vox = ceil(radius_mm ./ spacing);
    [Y, X, Z] = ndgrid(-radius_vox(1):radius_vox(1), ...
                        -radius_vox(2):radius_vox(2), ...
                        -radius_vox(3):radius_vox(3));
    dist_mm_sq = (Y*spacing(1)).^2 + (X*spacing(2)).^2 + (Z*spacing(3)).^2;
    se = strel(dist_mm_sq <= radius_mm^2);
    closed = imclose(mask, se);
end


function [marker_mask, artifact_w] = marker_and_artifact_maps(HU, marker_range, sigma_mm, spacing)
    lead = HU > 1200;
    flags = (HU >= marker_range(1) & HU <= marker_range(2)) & ...
            imdilate(lead, strel('sphere', 2));
    marker_mask = lead | flags;

    d_vox = bwdist(marker_mask);
    d_mm = d_vox * mean(spacing);
    artifact_w = exp(-(d_mm / sigma_mm).^2);
end


function bones = associate_tags(bones, tags)
    if isempty(tags) || isempty(bones), return; end
    centroids = zeros(numel(bones), 3);
    for i = 1:numel(bones)
        centroids(i,:) = bones{i}.centroid_mm;
    end

    fprintf('    Tag-bone distances (mm):\n');
    for t = 1:numel(tags)
        dists = vecnorm(centroids - tags{t}.centroid_mm, 2, 2);
        [d, idx] = min(dists);
        fprintf('      Marker %d -> Bone %d: %.1f mm\n', t, idx, d);
        if isempty(bones{idx}.tag_id) || d < bones{idx}.tag_dist
            bones{idx}.tag_id = tags{t}.label;
            bones{idx}.tag_dist = d;
        end
    end
end
