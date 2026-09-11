function out = run_bone_pipeline(dicomFolder, stlFolder, varargin)
% RUN_BONE_PIPELINE  Full pipeline for multi-bone CT segmentation and specimen packing.
%
%   out = run_bone_pipeline(dicomFolder, stlFolder)
%   out = run_bone_pipeline(dicomFolder, stlFolder, 'Option', value, ...)
%
% Pipeline stages:
%   1. DICOM loading
%   2. Bone separation (FMM envelope detection for excised-in-air specimens)
%   3. Cortical / cancellous segmentation (gradient-based)
%   4. Specimen packing (greedy mixed packing of STL shapes)
%   5. Visualization (3D overview)
%   6. Output saving (MAT + STL + NIfTI)
%
% Inputs
%   dicomFolder : path to folder containing DICOM CT series
%   stlFolder   : path to folder containing specimen STL files
%                  (Bend.STL, Compression.STL, Punch.STL, Shear.STL)
%                  May be '' when packing is disabled.
%
% Name-value options
%   'TagHUMin'            : 1200 (HU threshold for lead tag detection)
%   'MinBoneVolMM3'       : 500  (minimum bone component volume)
%   'ClosingRadiusMM'     : 3.0  (morphological closing radius)
%   'ArtifactSigmaMM'     : 3.0  (Gaussian falloff for artifact weighting)
%   'MaxBones'            : []   (keep only the N largest bones; [] = all)
%   'CorticalCancellous'  : true (run stage 3; false = bone mask only)
%   'PackSpecimens'       : true  (run specimen packing — slow)
%   'PackWholeBone'       : false (pack into full bone ignoring cortical/cancellous)
%   'PackingOrientations' : 6     (number of orientations per shape)
%   'PackingMinDepthMM'   : 0.5  (minimum depth for specimen placement)
%   'SaveOutputs'         : true (export MAT, NIfTI, STL files)
%   'SaveMat'             : true (include pipeline_results.mat — large)
%   'OutputDir'           : ''   (auto-create if empty)
%   'ShowViewer'          : true (show 3D visualization)
%
% Fallback for difficult scans (see bone.segment_options)
%   'Fallback'            : true  retry with a different preset when the
%                                 standard pass fails or looks contaminated
%   'FallbackPreset'      : ''    force a preset instead of choosing one
%                                 ('lowdensity' | 'tissue' | 'none')
%   'LowDensityHU'        : 250   below this interior median, flag as
%                                 osteoporotic
%   'MinContrastHU'       : 150   minimum bone-to-surroundings density step
%   'MaxLowFrac'          : 0.35  max fraction of mask below TissueCeilingHU
%   'TissueCeilingHU'     : 150   HU below which a voxel is not clearly bone
%   'MinFillRatio'        : 0.50  mask volume / source blob volume floor
%   'LowDensityFillRatio' : 0.75  retry a low-density bone below this fill
%
% Segmentation knobs (empty = preset default, see bone.segment_options)
%   'DensityScale', 'CorePrctile', 'FMMThreshMax', 'TissueScrub', 'MinBoneHU'
%
% Mesh export
%   'MeshSmoothIterations': 8    Taubin smoothing passes (0 = none)
%   'MeshSmoothLambda'    : 0.5  Taubin shrink weight
%   'MeshPreSmoothSigma'  : 1.0  Gaussian sigma applied to the mask first
%   'MeshDecimate'        : 0.5  fraction of faces kept in the smooth STL

% ---- Parse options ----
opts = struct( ...
    'TagHUMin',            1200, ...
    'MinBoneVolMM3',       500.0, ...
    'ClosingRadiusMM',     3.0, ...
    'ArtifactSigmaMM',     3.0, ...
    'MarkerRangeHU',       [200 700], ...
    'MaxBones',            [], ...
    'CorticalCancellous',  true, ...
    'PackSpecimens',       true, ...
    'PackWholeBone',       false, ...
    'PackingOrientations', 6, ...
    'PackingMinDepthMM',   0.5, ...
    'SaveOutputs',         true, ...
    'SaveMat',             true, ...
    'OutputDir',           '', ...
    'ShowViewer',          true, ...
    'TargetIsoMM',         [], ...
    'Smoothing',           false, ...
    'Fallback',            true, ...
    'FallbackPreset',      '', ...
    'LowDensityHU',        250, ...
    'MinContrastHU',       150, ...
    'MaxLowFrac',          0.35, ...
    'TissueCeilingHU',     150, ...
    'MinFillRatio',        0.50, ...
    'LowDensityFillRatio', 0.75, ...
    'DensityScale',        [], ...
    'CorePrctile',         [], ...
    'FMMThreshMax',        [], ...
    'TissueScrub',         [], ...
    'MinBoneHU',           [], ...
    'MeshSmoothIterations', 8, ...
    'MeshSmoothLambda',    0.5, ...
    'MeshPreSmoothSigma',  1.0, ...
    'MeshDecimate',        0.5 ...
);
opts = utils.parse_opts(opts, varargin{:});

if nargin < 2, stlFolder = ''; end

% Cortical/cancellous masks are the input to packing — without them there is
% nothing to pack into.
if ~opts.CorticalCancellous && opts.PackSpecimens
    opts.PackSpecimens = false;
end

t_start = tic;

fprintf('\n');
fprintf('==========================================================\n');
fprintf('  BONE SEGMENTATION PIPELINE\n');
fprintf('==========================================================\n');
fprintf('  DICOM : %s\n', dicomFolder);
if isempty(stlFolder)
    fprintf('  STL   : (none)\n');
else
    fprintf('  STL   : %s\n', stlFolder);
end
if ~opts.CorticalCancellous
    fprintf('  Mode  : segmentation only (no cortical/cancellous, no packing)\n');
end
fprintf('==========================================================\n\n');

% ==== Stage 1: DICOM Loading ====
fprintf('[1/6] Loading DICOM series...\n');
t1 = tic;
ds = dicom.series_load(dicomFolder, ...
    'TargetIsoMM', opts.TargetIsoMM, 'Smoothing', opts.Smoothing);
fprintf('[1/6] Done (%.1fs)\n', toc(t1));
fprintf('       Volume %dx%dx%d  |  spacing [%.3f  %.3f  %.3f] mm  |  HU [%.0f, %.0f]\n\n', ...
    ds.size(1), ds.size(2), ds.size(3), ds.spacing, min(ds.HU(:)), max(ds.HU(:)));

% ==== Stage 2: Bone Separation ====
fprintf('[2/6] Separating bones...');
t2 = tic;
sep_result = bone.separate_bones(ds, bone.segment_options(opts));
sep_result.preset = 'standard';
sep_result = attach_quality(sep_result, ds, opts);
fprintf(' done (%.1fs)\n', toc(t2));
fprintf('       %d bones found, %d markers detected\n', ...
    numel(sep_result.bones), sep_result.n_tags);
print_quality(sep_result);

% ---- Fallback pass for low-density or tissue-contaminated scans ----
% The standard preset is tuned for healthy cortical bone. When it comes back
% empty, short, or with a mask that does not look separable from what
% surrounds it, retry under different assumptions and keep the better of the
% two. Ties go to the standard pass.
preset = fallback_preset(sep_result, opts);
if ~isempty(preset)
    fprintf('       Fallback: retrying with ''%s'' preset\n', preset);
    t2b = tic;
    alt_opts = bone.segment_options(opts, preset);
    try
        sep_alt = bone.separate_bones(ds, alt_opts);
        sep_alt.preset = preset;
        sep_alt = attach_quality(sep_alt, ds, opts);
        [take, why] = accept_fallback(sep_result, sep_alt, preset);
        fprintf('       Fallback (%.1fs): %s — %s\n', toc(t2b), ...
            ternary(take, 'ACCEPTED', 'rejected'), why);
        if take
            sep_result = sep_alt;
            print_quality(sep_result);
        end
    catch ME
        fprintf('       Fallback failed (%s) — keeping standard result\n', ME.message);
    end
end
fprintf('\n');

n_bones = numel(sep_result.bones);

if n_bones == 0
    warning('No bones found. Check DICOM data and thresholds.');
    out = struct('ds', ds, 'separation', sep_result, ...
        'segmentation', {{}}, 'packing', {{}});
    return;
end

% Keep only the N largest bones (separate_bones sorts by volume descending).
% For single-bone scans set MaxBones = 1 to drop spurious extra objects.
% n_bones_found records what separation actually produced, so a scan where
% extras were discarded is still visible after the fact.
sep_result.n_bones_found = n_bones;
if ~isempty(opts.MaxBones) && n_bones > opts.MaxBones
    fprintf('       Keeping %d largest of %d bones (MaxBones)\n\n', opts.MaxBones, n_bones);
    sep_result.bones = sep_result.bones(1:opts.MaxBones);
    n_bones = numel(sep_result.bones);
end

bone_masks = cell(1, n_bones);
for bi = 1:n_bones
    bone_masks{bi} = sep_result.bones{bi}.mask;
end

use_parallel = ~isempty(ver('parallel')) && n_bones > 1;

% ==== Stage 3: Cortical / Cancellous Segmentation ====
seg_results = cell(1, n_bones);
if ~opts.CorticalCancellous
    fprintf('[3/6] Cortical/cancellous segmentation skipped (CorticalCancellous=false)\n\n');
else
    fprintf('[3/6] Cortical/cancellous segmentation...');
    t3 = tic;
    if use_parallel
        parfor bi = 1:n_bones
            [cort, canc, seg_info] = bone.cortical_cancellous(ds, bone_masks{bi}, opts);
            seg_results{bi} = struct('cortical', cort, 'cancellous', canc, 'info', seg_info);
        end
    else
        for bi = 1:n_bones
            [cort, canc, seg_info] = bone.cortical_cancellous(ds, bone_masks{bi}, opts);
            seg_results{bi} = struct('cortical', cort, 'cancellous', canc, 'info', seg_info);
        end
    end
    fprintf(' done (%.1fs)%s\n\n', toc(t3), ternary(use_parallel, ' [parallel]', ''));
end
has_seg = opts.CorticalCancellous;

% ==== Stage 4: Specimen Packing ====
pack_results = cell(1, n_bones);
if ~opts.PackSpecimens
    fprintf('[4/6] Specimen packing skipped (PackSpecimens=false)\n\n');
else
    stl_names = {'Bend', 'Compression', 'Punch', 'Shear'};
    stl_paths = {};
    stl_found = {};
    for si = 1:numel(stl_names)
        candidates = {fullfile(stlFolder, [stl_names{si} '.STL']), ...
                      fullfile(stlFolder, [stl_names{si} '.stl'])};
        for ci = 1:numel(candidates)
            if exist(candidates{ci}, 'file')
                stl_paths{end+1} = candidates{ci}; %#ok<AGROW>
                stl_found{end+1} = stl_names{si}; %#ok<AGROW>
                break;
            end
        end
    end

    if isempty(stl_paths)
        fprintf('[4/6] Specimen packing skipped (no STL files found)\n\n');
    else
        if opts.PackWholeBone
            fprintf('[4/6] Packing specimens — whole bone (%s)\n', strjoin(stl_found, ', '));
        else
            fprintf('[4/6] Packing specimens — cortical/cancellous (%s)\n', strjoin(stl_found, ', '));
        end
        t4 = tic;

        % Print specimen dimensions (load once for display)
        fprintf('       Specimen dimensions:\n');
        fprintf('       %-15s  %6s x %6s x %6s  %10s\n', 'Shape', 'X mm', 'Y mm', 'Z mm', 'Volume');
        fprintf('       %-15s  %6s   %6s   %6s  %10s\n', '---------------', '------', '------', '------', '----------');
        for si = 1:numel(stl_paths)
            try
                TR = stlread(stl_paths{si});
                V_raw = double(TR.Points);
                V_raw = V_raw - mean(V_raw, 1);
                bbox = max(V_raw, [], 1) - min(V_raw, [], 1);
                % Approximate volume from convex hull
                try
                    [~, stl_vol] = convhull(V_raw);
                catch
                    stl_vol = prod(bbox);
                end
                fprintf('       %-15s  %5.1f  x %5.1f  x %5.1f   %7.0f mm3\n', ...
                    stl_found{si}, bbox(1), bbox(2), bbox(3), stl_vol);
            catch
                fprintf('       %-15s  (failed to read)\n', stl_found{si});
            end
        end

        % Print bone region volumes for comparison
        fprintf('       Bone volumes:\n');
        for bi = 1:n_bones
            bone_vol_bi = sep_result.bones{bi}.volume_mm3;
            if opts.PackWholeBone
                fprintf('       Bone #%d: %.0f mm3 (whole)\n', bi, bone_vol_bi);
            else
                cort_v = seg_results{bi}.info.cortical_volume_mm3;
                canc_v = seg_results{bi}.info.cancellous_volume_mm3;
                fprintf('       Bone #%d: %.0f mm3 (cortical %.0f, cancellous %.0f)\n', ...
                    bi, bone_vol_bi, cort_v, canc_v);
            end
        end
        fprintf('\n');

        bone_axes = cell(1, n_bones);
        corticals = cell(1, n_bones);
        cancellouses = cell(1, n_bones);
        for bi = 1:n_bones
            bone_axes{bi} = [0; 0; 1];
            if isfield(seg_results{bi}.info, 'bone_shape')
                bm = bone_masks{bi};
                [rr, cc, ss] = ind2sub(size(ds.HU), find(bm));
                coords = [rr(:)*ds.spacing(1), cc(:)*ds.spacing(2), ss(:)*ds.spacing(3)];
                coords = coords - mean(coords, 1);
                [V, ~] = eig(coords' * coords);
                bone_axes{bi} = V(:, 3);
            end
            corticals{bi} = seg_results{bi}.cortical;
            cancellouses{bi} = seg_results{bi}.cancellous;
        end

        if use_parallel
            parfor bi = 1:n_bones
                pack_results{bi} = bone.pack_specimens( ...
                    bone_masks{bi}, corticals{bi}, cancellouses{bi}, ...
                    ds, stl_paths, stl_found, opts, bone_axes{bi}, bi);
            end
        else
            for bi = 1:n_bones
                pack_results{bi} = bone.pack_specimens( ...
                    bone_masks{bi}, corticals{bi}, cancellouses{bi}, ...
                    ds, stl_paths, stl_found, opts, bone_axes{bi}, bi);
            end
        end
        fprintf('       Done (%.1fs)%s\n\n', toc(t4), ternary(use_parallel, ' [parallel]', ''));
    end
end

% ==== Stage 5: Visualization ====
if opts.ShowViewer && ~has_seg
    fprintf('[5/6] Visualization skipped (needs cortical/cancellous segmentation)\n\n');
elseif opts.ShowViewer
    fprintf('[5/6] Generating visualizations...');
    t5 = tic;
    bone.visualize_results(ds, sep_result, seg_results, pack_results, opts);
    fprintf(' done (%.1fs)\n\n', toc(t5));
else
    fprintf('[5/6] Visualization skipped\n\n');
end

% ==== Stage 6: Save Outputs ====
if opts.SaveOutputs
    fprintf('[6/6] Saving outputs...');
    t6 = tic;

    if isempty(opts.OutputDir)
        [parentDir, seriesName] = fileparts(string(dicomFolder));
        baseOut = fullfile(parentDir, 'bone_pipeline_outputs', seriesName);
        tstamp = datestr(now, 'yyyymmdd_HHMMSS');
        outDir = fullfile(baseOut, tstamp);
    else
        outDir = opts.OutputDir;
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    opts.OutputDir = outDir;

    % --- MAT file ---
    if opts.SaveMat
        try
            save(fullfile(outDir, 'pipeline_results.mat'), ...
                'sep_result', 'seg_results', 'pack_results', 'opts', '-v7.3');
        catch ME
            warning('Save MAT failed: %s', ME.message);
        end
    end

    % --- Per-bone NIfTI + STL ---
    mesh_vols = zeros(1, n_bones);
    for bi = 1:n_bones
        bm = sep_result.bones{bi}.mask;

        % NIfTI masks
        try
            write_mask_nifti(fullfile(outDir, sprintf('bone_%02d_mask.nii.gz', bi)), bm, ds);
            if has_seg
                write_mask_nifti(fullfile(outDir, sprintf('bone_%02d_cortical.nii.gz', bi)), ...
                    seg_results{bi}.cortical, ds);
                write_mask_nifti(fullfile(outDir, sprintf('bone_%02d_cancellous.nii.gz', bi)), ...
                    seg_results{bi}.cancellous, ds);
            end

            HU_masked = int16(ds.HU);
            HU_masked(~bm) = -3000;
            write_volume_nifti(fullfile(outDir, sprintf('bone_%02d_hu.nii.gz', bi)), HU_masked, ds);
        catch ME
            warning('NIfTI save failed for bone %d: %s', bi, ME.message);
        end

        % Voxelized STL (fast, preserves exact mask geometry)
        try
            if any(bm(:))
                fv = isosurface(smooth3(double(bm), 'gaussian', 3), 0.5);
                if ~isempty(fv.vertices)
                    fv.vertices(:,1) = fv.vertices(:,1) * ds.spacing(2);
                    fv.vertices(:,2) = fv.vertices(:,2) * ds.spacing(1);
                    fv.vertices(:,3) = fv.vertices(:,3) * ds.spacing(3);
                    mesh = struct('vertices', fv.vertices, 'faces', fv.faces);
                    meshing.write_stl_binary(fullfile(outDir, sprintf('bone_%02d_voxelized.stl', bi)), mesh);
                end
            end
        catch ME
            warning('Voxelized STL save failed for bone %d: %s', bi, ME.message);
        end

        % Smooth anatomical STL (Taubin-smoothed, decimated).
        % Taubin alternates a shrinking and an inflating pass, so the
        % surface loses its voxel staircase without pulling in off the
        % bone the way repeated Laplacian passes do.
        try
            if any(bm(:))
                smooth_field = smooth3(double(bm), 'gaussian', 5, opts.MeshPreSmoothSigma);
                fv_smooth = isosurface(smooth_field, 0.5);
                if ~isempty(fv_smooth.vertices) && size(fv_smooth.faces, 1) > 100
                    fv_smooth.vertices(:,1) = fv_smooth.vertices(:,1) * ds.spacing(2);
                    fv_smooth.vertices(:,2) = fv_smooth.vertices(:,2) * ds.spacing(1);
                    fv_smooth.vertices(:,3) = fv_smooth.vertices(:,3) * ds.spacing(3);

                    V_s = meshing.smooth_mesh_taubin(fv_smooth.vertices, ...
                        fv_smooth.faces, opts.MeshSmoothIterations, ...
                        opts.MeshSmoothLambda);

                    smooth_mesh = struct('vertices', V_s, 'faces', fv_smooth.faces);
                    meshing.write_stl_binary( ...
                        fullfile(outDir, sprintf('bone_%02d_smooth.stl', bi)), ...
                        smooth_mesh, 'Decimate', opts.MeshDecimate);

                    % How far the smoothing moved the surface, as a check
                    % that it is not quietly eating the bone.
                    mesh_vols(bi) = mesh_volume(V_s, fv_smooth.faces);
                end
            end
        catch ME
            warning('Smooth STL save failed for bone %d: %s', bi, ME.message);
        end
    end


    % --- Text summary ---
    try
        write_summary_file(fullfile(outDir, 'pipeline_summary.txt'), ...
            ds, sep_result, seg_results, pack_results, dicomFolder, n_bones, has_seg);
    catch ME
        warning('Summary save failed: %s', ME.message);
    end

    fprintf(' done (%.1fs)\n', toc(t6));

    % How far smoothing moved the surface, as a check that it is not
    % quietly eating the bone.
    for bi = 1:n_bones
        if mesh_vols(bi) > 0
            mask_vol = sep_result.bones{bi}.volume_mm3;
            fprintf('       Bone #%d smooth STL %.0f mm3 vs mask %.0f mm3 (%+.1f%%)\n', ...
                bi, mesh_vols(bi), mask_vol, ...
                (mesh_vols(bi) / max(mask_vol, eps) - 1) * 100);
        end
    end

    fprintf('       Output: %s\n\n', outDir);
    out.outputDir = outDir;
    out.meshVolumes = mesh_vols;
else
    fprintf('[6/6] Output saving skipped\n\n');
end

% ==== Build output struct ====
out.ds = ds;
out.separation = sep_result;
out.segmentation = seg_results;
out.packing = pack_results;

% ==== Final Summary Table ====
elapsed = toc(t_start);

fprintf('==========================================================\n');
fprintf('  PIPELINE COMPLETE  (%.1f s)\n', elapsed);
fprintf('==========================================================\n\n');

% --- Bone summary table ---
total_vol = 0;
total_cort = 0;
total_canc = 0;

if has_seg
    fprintf('  %-6s  %10s  %8s  %8s  %10s  %10s  %8s  %5s\n', ...
        'Bone', 'Volume', 'Mean HU', 'Shape', 'Cortical', 'Cancellous', 'Cort %%', 'Tag');
    fprintf('  %-6s  %10s  %8s  %8s  %10s  %10s  %8s  %5s\n', ...
        '------', '----------', '--------', '--------', '----------', '----------', '--------', '-----');
else
    fprintf('  %-6s  %10s  %8s  %5s\n', 'Bone', 'Volume', 'Mean HU', 'Tag');
    fprintf('  %-6s  %10s  %8s  %5s\n', '------', '----------', '--------', '-----');
end

for bi = 1:n_bones
    b = sep_result.bones{bi};
    total_vol = total_vol + b.volume_mm3;

    if ~isempty(b.tag_id)
        tag_str = sprintf('%d', b.tag_id);
    else
        tag_str = '-';
    end

    if has_seg
        si = seg_results{bi}.info;
        total_cort = total_cort + si.cortical_volume_mm3;
        total_canc = total_canc + si.cancellous_volume_mm3;
        fprintf('  %-6s  %8.0f mm3  %6.0f HU  %8s  %8.0f mm3  %8.0f mm3  %6.1f%%  %5s\n', ...
            sprintf('#%d', bi), b.volume_mm3, b.mean_hu, si.bone_shape, ...
            si.cortical_volume_mm3, si.cancellous_volume_mm3, ...
            si.cortical_fraction * 100, tag_str);
    else
        fprintf('  %-6s  %8.0f mm3  %6.0f HU  %5s\n', ...
            sprintf('#%d', bi), b.volume_mm3, b.mean_hu, tag_str);
    end
end

if has_seg
    fprintf('  %-6s  %8.0f mm3  %8s  %8s  %8.0f mm3  %8.0f mm3  %6.1f%%\n', ...
        'TOTAL', total_vol, '', '', total_cort, total_canc, ...
        total_cort / max(1, total_cort + total_canc) * 100);
else
    fprintf('  %-6s  %8.0f mm3\n', 'TOTAL', total_vol);
end

% --- Segmentation quality ---
fprintf('\n  Segmentation pass: %s\n', sep_result.preset);
needs_review = false;
for bi = 1:n_bones
    if ~isfield(sep_result.bones{bi}, 'quality'), continue; end
    q = sep_result.bones{bi}.quality;
    if isempty(q.flags)
        flag_str = 'ok';
    else
        flag_str = strjoin(q.flags, ', ');
        needs_review = true;
    end
    fprintf('  #%d  core %.0f HU  rind %.0f HU  contrast %.0f  %.0f%% low  ->  %s\n', ...
        bi, q.core_hu, q.rind_hu, q.contrast, q.low_frac * 100, flag_str);
end
if needs_review
    fprintf('\n  Flagged bones are the best this data supports — check them by eye\n');
    fprintf('  before using them; the flag says why, not that the mask is wrong.\n');
end

% --- Packing summary table ---
has_packing = ~isempty(pack_results) && ~isempty(pack_results{1}) && ...
    isstruct(pack_results{1}) && isfield(pack_results{1}, 'n_total');
if has_packing
    is_whole = isfield(pack_results{1}, 'whole_bone') && pack_results{1}.whole_bone;

    if is_whole
        fprintf('\n  Packing mode: whole bone\n');
        fprintf('  %-6s  %8s  %8s\n', 'Bone', 'Placed', 'Total');
        fprintf('  %-6s  %8s  %8s\n', '------', '--------', '--------');

        grand_total = 0;
        for bi = 1:n_bones
            pr = pack_results{bi};
            grand_total = grand_total + pr.n_total;
            fprintf('  %-6s  %8d  %8d\n', sprintf('#%d', bi), pr.n_whole, pr.n_total);
        end
        fprintf('  %-6s  %8d  %8d\n', 'TOTAL', grand_total, grand_total);
    else
        fprintf('\n  Packing mode: cortical / cancellous\n');
        fprintf('  %-6s  %8s  %8s  %8s\n', 'Bone', 'Cortical', 'Cancel.', 'Total');
        fprintf('  %-6s  %8s  %8s  %8s\n', '------', '--------', '--------', '--------');

        grand_cort = 0; grand_canc = 0;
        for bi = 1:n_bones
            pr = pack_results{bi};
            grand_cort = grand_cort + pr.n_cortical;
            grand_canc = grand_canc + pr.n_cancellous;
            fprintf('  %-6s  %8d  %8d  %8d\n', sprintf('#%d', bi), ...
                pr.n_cortical, pr.n_cancellous, pr.n_total);
        end
        fprintf('  %-6s  %8d  %8d  %8d\n', 'TOTAL', ...
            grand_cort, grand_canc, grand_cort + grand_canc);
    end

    % Per-type breakdown
    if isfield(pack_results{1}, 'summary')
        fprintf('\n  Specimens by type:\n');
        all_types = fieldnames(pack_results{1}.summary);
        for ti = 1:numel(all_types)
            type_total = 0;
            for bi = 1:n_bones
                if isfield(pack_results{bi}.summary, all_types{ti})
                    s = pack_results{bi}.summary.(all_types{ti});
                    type_total = type_total + s.cortical + s.cancellous + s.whole;
                end
            end
            fprintf('    %-15s : %d\n', all_types{ti}, type_total);
        end
    end
end

fprintf('\n==========================================================\n\n');
end


% =========================================================================
%  Local helper functions
% =========================================================================

function sep = attach_quality(sep, ds, opts)
% Score every bone against the material around it.
    for bi = 1:numel(sep.bones)
        qo = struct( ...
            'LowDensityHU',    opts.LowDensityHU, ...
            'MinContrastHU',   opts.MinContrastHU, ...
            'MaxLowFrac',      opts.MaxLowFrac, ...
            'TissueCeilingHU', opts.TissueCeilingHU, ...
            'MinFillRatio',    opts.MinFillRatio, ...
            'SourceVolMM3',    0);
        if isfield(sep.bones{bi}, 'source_vol_mm3')
            qo.SourceVolMM3 = sep.bones{bi}.source_vol_mm3;
        end
        sep.bones{bi}.quality = bone.mask_quality( ...
            ds.HU, sep.bones{bi}.mask, ds.spacing, qo);
    end
end


function print_quality(sep)
    for bi = 1:numel(sep.bones)
        q = sep.bones{bi}.quality;
        if isempty(q.flags)
            flag_str = 'ok';
        else
            flag_str = strjoin(q.flags, ', ');
        end
        fprintf('       Bone #%d [%s]: core %.0f HU, rind %.0f HU, contrast %.0f, %.0f%% low — %s\n', ...
            bi, sep.preset, q.core_hu, q.rind_hu, q.contrast, q.low_frac * 100, flag_str);
    end
end


function preset = fallback_preset(sep, opts)
% Decide whether a second pass is worth running, and under which assumption.
% Low density on its own is not a failure — plenty of osteoporotic bone
% segments cleanly — so it only triggers a retry when the mask also came
% back noticeably short of the blob it grew from.
    preset = '';
    if ~opts.Fallback, return; end

    if ~isempty(opts.FallbackPreset)
        if strcmpi(opts.FallbackPreset, 'none'), return; end
        preset = lower(char(opts.FallbackPreset));
        return;
    end

    if isempty(sep.bones)
        preset = 'lowdensity';
        return;
    end

    q = sep.bones{1}.quality;   % the largest bone decides
    if any(strcmp(q.flags, 'tissue_suspect'))
        preset = 'tissue';
    elseif any(strcmp(q.flags, 'under_segmented'))
        preset = 'lowdensity';
    elseif any(strcmp(q.flags, 'low_density')) && ...
           q.fill_ratio > 0 && q.fill_ratio < opts.LowDensityFillRatio
        preset = 'lowdensity';
    end
end


function [take, why] = accept_fallback(sep_std, sep_alt, preset)
% Keep the fallback only when it clearly beats the standard pass. Anything
% ambiguous keeps the standard result — the point of the fallback is to
% rescue scans that failed, not to second-guess ones that worked.
    take = false;

    if isempty(sep_alt.bones)
        why = 'fallback found no bone';
        return;
    end
    if isempty(sep_std.bones)
        take = true;
        why = 'standard pass found no bone';
        return;
    end

    a = sep_alt.bones{1};  qa = a.quality;
    s = sep_std.bones{1};  qs = s.quality;
    vol_ratio = a.volume_mm3 / max(s.volume_mm3, eps);

    switch preset
        case 'lowdensity'
            grew = vol_ratio > 1.05;
            % A bigger mask that is no longer separable from its
            % surroundings is tissue, not recovered bone.
            held_up = qa.score >= qs.score - 0.05;
            new_tissue = any(strcmp(qa.flags, 'tissue_suspect')) && ...
                        ~any(strcmp(qs.flags, 'tissue_suspect'));
            take = grew && held_up && ~new_tissue;
            if new_tissue
                why = sprintf('%.0f%% larger but now tissue-suspect', (vol_ratio-1)*100);
            elseif ~grew
                why = sprintf('no larger than standard (%.0f%%)', vol_ratio*100);
            elseif ~held_up
                why = sprintf('grew %.0f%% but score fell %.2f to %.2f', ...
                    (vol_ratio-1)*100, qs.score - qa.score, qa.score);
            else
                why = sprintf('recovered %.0f%% more bone, score %.2f vs %.2f', ...
                    (vol_ratio-1)*100, qa.score, qs.score);
            end

        case 'tissue'
            cleaner = qa.score > qs.score + 0.05;
            % Stripping tissue should trim the mask, not gut it.
            intact = vol_ratio > 0.5;
            take = cleaner && intact;
            if ~cleaner
                why = sprintf('no cleaner (score %.2f vs %.2f)', qa.score, qs.score);
            elseif ~intact
                why = sprintf('removed too much (%.0f%% of standard)', vol_ratio*100);
            else
                why = sprintf('score %.2f vs %.2f, kept %.0f%% of volume', ...
                    qa.score, qs.score, vol_ratio*100);
            end

        otherwise
            take = qa.score > qs.score + 0.05;
            why = sprintf('score %.2f vs %.2f', qa.score, qs.score);
    end
end


function v = mesh_volume(V, F)
% Enclosed volume of a closed triangle mesh (divergence theorem).
    v = 0;
    if isempty(V) || isempty(F), return; end
    v1 = V(F(:,1), :);  v2 = V(F(:,2), :);  v3 = V(F(:,3), :);
    v = abs(sum(dot(v1, cross(v2, v3, 2), 2))) / 6;
end


function write_mask_nifti(filename, mask, ds)
    data = int16(mask) * 1000;
    write_volume_nifti(filename, data, ds);
end


function write_volume_nifti(filename, data, ds)
    M = [ds.dir_row'   * ds.spacing(1), ...
         ds.dir_col'   * ds.spacing(2), ...
         ds.dir_slice' * ds.spacing(3), ...
         ds.origin'; ...
         0 0 0 1];

    seedFile = fullfile(tempdir, 'bone_pipeline_seed.nii');
    if exist(seedFile, 'file'), delete(seedFile); end
    niftiwrite(zeros(size(data), 'int16'), seedFile);
    info = niftiinfo(seedFile);
    delete(seedFile);

    info.Datatype = 'int16';
    info.PixelDimensions = ds.spacing;
    info.Transform = affine3d(M');

    niftiwrite(int16(data), filename, info, 'Compressed', true);
end


function write_summary_file(filepath, ds, sep_result, seg_results, pack_results, dicomFolder, n_bones, has_seg)
    fid = fopen(filepath, 'w');
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, 'BONE SEGMENTATION PIPELINE SUMMARY\n');
    fprintf(fid, '===================================\n');
    fprintf(fid, 'Date   : %s\n', datestr(now));
    fprintf(fid, 'DICOM  : %s\n', dicomFolder);
    fprintf(fid, 'Volume : %dx%dx%d, spacing [%.3f %.3f %.3f] mm\n', ...
        ds.size(1), ds.size(2), ds.size(3), ds.spacing);
    fprintf(fid, 'HU range: [%.0f, %.0f]\n\n', min(ds.HU(:)), max(ds.HU(:)));
    fprintf(fid, 'Bones   : %d\n', n_bones);
    fprintf(fid, 'Markers : %d\n', sep_result.n_tags);
    if isfield(sep_result, 'preset')
        fprintf(fid, 'Pass    : %s\n', sep_result.preset);
    end
    fprintf(fid, '\n');

    for bi = 1:n_bones
        b = sep_result.bones{bi};
        if ~isempty(b.tag_id)
            tag_str = sprintf('tag %d (%.1f mm)', b.tag_id, b.tag_dist);
        else
            tag_str = 'no tag';
        end
        if has_seg
            si = seg_results{bi}.info;
            fprintf(fid, 'Bone %d: %.1f mm3 | HU %.0f | %s | %s\n', ...
                bi, b.volume_mm3, b.mean_hu, si.bone_shape, tag_str);
            fprintf(fid, '  Cortical: %.0f mm3 (%.1f%%) depth %.2f mm\n', ...
                si.cortical_volume_mm3, si.cortical_fraction*100, si.mean_cortical_depth_mm);
            fprintf(fid, '  Cancellous: %.0f mm3\n', si.cancellous_volume_mm3);
        else
            fprintf(fid, 'Bone %d: %.1f mm3 | HU %.0f | %s\n', ...
                bi, b.volume_mm3, b.mean_hu, tag_str);
        end
        if isfield(b, 'quality')
            q = b.quality;
            if isempty(q.flags), flag_str = 'ok'; else, flag_str = strjoin(q.flags, ', '); end
            fprintf(fid, '  Quality: core %.0f HU, rind %.0f HU, contrast %.0f, %.0f%% low-HU -> %s\n', ...
                q.core_hu, q.rind_hu, q.contrast, q.low_frac * 100, flag_str);
        end
    end

    has_packing = ~isempty(pack_results) && ~isempty(pack_results{1}) && ...
        isstruct(pack_results{1}) && isfield(pack_results{1}, 'n_total');
    if has_packing
        is_whole = isfield(pack_results{1}, 'whole_bone') && pack_results{1}.whole_bone;
        fprintf(fid, '\nSPECIMEN PACKING (%s)\n', ternary(is_whole, 'whole bone', 'cortical/cancellous'));
        fprintf(fid, '================\n');
        for bi = 1:n_bones
            pr = pack_results{bi};
            if is_whole
                fprintf(fid, 'Bone %d: %d specimens (whole bone)\n', bi, pr.n_total);
            else
                fprintf(fid, 'Bone %d: %d cortical + %d cancellous = %d specimens\n', ...
                    bi, pr.n_cortical, pr.n_cancellous, pr.n_total);
            end
            if isfield(pr, 'summary')
                fnames = fieldnames(pr.summary);
                for fi = 1:numel(fnames)
                    s = pr.summary.(fnames{fi});
                    if is_whole
                        fprintf(fid, '  %-15s: %d\n', fnames{fi}, s.whole);
                    else
                        fprintf(fid, '  %-15s: %d cortical, %d cancellous\n', ...
                            fnames{fi}, s.cortical, s.cancellous);
                    end
                end
            end
        end
    end
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
