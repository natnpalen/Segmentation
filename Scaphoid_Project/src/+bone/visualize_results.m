function visualize_results(ds, sep_result, seg_results, pack_results, opts)
% VISUALIZE_RESULTS  Diagnostic visualization for the bone pipeline.
%
%   bone.visualize_results(ds, sep_result, seg_results, pack_results, opts)
%
% Creates:
%   Figure 1: 3D overview — all bones color-coded with markers
%   Figure 2: Axial slice montage — bone masks overlaid on CT
%   Figure 3: Per-bone HU histograms (marker contamination check)

vol = double(ds.HU);
spacing = ds.spacing;
bones = sep_result.bones;
n_bones = numel(bones);

colors = [0.2 0.6 1.0;   % blue
          1.0 0.5 0.0;   % orange
          0.2 0.8 0.3;   % green
          0.9 0.2 0.9;   % magenta
          0.0 0.8 0.8;   % cyan
          1.0 0.8 0.2;   % yellow
          0.6 0.3 0.1;   % brown
          0.5 0.5 0.5];  % gray
if n_bones > size(colors, 1)
    colors = [colors; lines(n_bones - size(colors, 1))];
end

% ========================================================================
%  Figure 1: 3D overview — all bones + markers
% ========================================================================
fig1 = figure('Name', 'Bone Separation Overview', 'Color', 'w', ...
    'Position', [50 50 900 700]);

for bi = 1:n_bones
    mask_i = bones{bi}.mask;
    if ~any(mask_i(:)), continue; end

    try
        fv = isosurface(smooth3(double(mask_i), 'gaussian', 3), 0.5);
        if isempty(fv.vertices), continue; end
        fv.vertices(:,1) = fv.vertices(:,1) * spacing(2);
        fv.vertices(:,2) = fv.vertices(:,2) * spacing(1);
        fv.vertices(:,3) = fv.vertices(:,3) * spacing(3);

        patch(fv, 'FaceColor', colors(bi,:), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.7);
        hold on;

        cm = bones{bi}.centroid_mm;
        if ~isempty(bones{bi}.tag_id)
            lbl = sprintf('Bone %d (tag %d)\n%.0f mm^3, HU=%.0f', ...
                bi, bones{bi}.tag_id, bones{bi}.volume_mm3, bones{bi}.mean_hu);
        else
            lbl = sprintf('Bone %d\n%.0f mm^3, HU=%.0f', ...
                bi, bones{bi}.volume_mm3, bones{bi}.mean_hu);
        end
        text(cm(2), cm(1), cm(3), lbl, 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', colors(bi,:)*0.6);
    catch
        continue;
    end
end

% Show markers as red
if isfield(sep_result, 'marker_mask') && any(sep_result.marker_mask(:))
    mk = sep_result.marker_mask;
    try
        fv_mk = isosurface(smooth3(double(mk), 'gaussian', 3), 0.5);
        if ~isempty(fv_mk.vertices)
            fv_mk.vertices(:,1) = fv_mk.vertices(:,1) * spacing(2);
            fv_mk.vertices(:,2) = fv_mk.vertices(:,2) * spacing(1);
            fv_mk.vertices(:,3) = fv_mk.vertices(:,3) * spacing(3);
            patch(fv_mk, 'FaceColor', [1 0 0], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.9);
        end
    catch
    end
end

axis equal vis3d off;
camlight headlight; lighting gouraud;
title(sprintf('Bone Separation: %d bones found', n_bones));
rotate3d on;

% ========================================================================
%  Figure 2: Axial slice montage — bone masks overlaid on CT
% ========================================================================
fig2 = figure('Name', 'Bone Separation — Axial Slices', 'Color', 'k', ...
    'Position', [100 100 1400 800]);

% Pick slices that pass through bone tissue
bone_slices = [];
for bi = 1:n_bones
    [~, ~, ss] = ind2sub(size(vol), find(bones{bi}.mask));
    bone_slices = [bone_slices; unique(ss)]; %#ok<AGROW>
end
bone_slices = unique(bone_slices);
if isempty(bone_slices)
    bone_slices = round(linspace(1, size(vol, 3), 16));
end

% Sample ~16 evenly spaced slices from the bone range
n_slices = min(16, numel(bone_slices));
idx = round(linspace(1, numel(bone_slices), n_slices));
show_slices = bone_slices(idx);

n_rows = ceil(n_slices / 4);
n_cols = min(4, n_slices);

for si = 1:n_slices
    subplot(n_rows, n_cols, si);

    slice_idx = show_slices(si);
    ct_slice = vol(:, :, slice_idx);

    % Normalize CT to [0,1] for display (window: -200 to 1500 HU)
    ct_disp = (ct_slice - (-200)) / (1500 - (-200));
    ct_disp = max(0, min(1, ct_disp));

    % Create RGB image from grayscale CT
    rgb = repmat(ct_disp, [1 1 3]);

    % Overlay each bone mask in its color
    for bi = 1:n_bones
        bone_slice = bones{bi}.mask(:, :, slice_idx);
        if ~any(bone_slice(:)), continue; end

        perim = bwperim(bone_slice);
        thick_perim = imdilate(perim, strel('disk', 1));

        for ch = 1:3
            layer = rgb(:,:,ch);
            % Fill interior with tinted overlay
            layer(bone_slice) = layer(bone_slice) * 0.6 + colors(bi, ch) * 0.4;
            % Bright boundary
            layer(thick_perim) = colors(bi, ch);
            rgb(:,:,ch) = layer;
        end
    end

    % Overlay marker mask in red
    if isfield(sep_result, 'marker_mask')
        mk_slice = sep_result.marker_mask(:, :, slice_idx);
        if any(mk_slice(:))
            mk_perim = imdilate(bwperim(mk_slice), strel('disk', 1));
            rgb(:,:,1) = max(rgb(:,:,1), double(mk_perim));
            rgb(:,:,2) = rgb(:,:,2) .* (1 - double(mk_perim)*0.8);
            rgb(:,:,3) = rgb(:,:,3) .* (1 - double(mk_perim)*0.8);
        end
    end

    imshow(rgb);
    title(sprintf('z=%d (%.1f mm)', slice_idx, slice_idx * spacing(3)), ...
        'Color', 'w', 'FontSize', 8);
end

sgtitle('Bone Masks on CT — Axial Slices', 'Color', 'w', 'FontSize', 12);

% ========================================================================
%  Figure 3: Per-bone HU histograms
% ========================================================================
fig3 = figure('Name', 'Per-Bone HU Distribution', 'Color', 'w', ...
    'Position', [150 150 1200 400]);

n_cols_hist = min(n_bones, 5);
for bi = 1:min(n_bones, n_cols_hist)
    subplot(1, n_cols_hist, bi);

    hu_vals = vol(bones{bi}.mask);

    histogram(hu_vals, 80, 'FaceColor', colors(bi,:), 'EdgeColor', 'none', ...
        'FaceAlpha', 0.8);
    hold on;

    % Mark key thresholds
    xline(0, '--', 'Air/Tissue', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xline(1200, '--', 'Metal', 'Color', [1 0 0], 'LineWidth', 1);
    if bones{bi}.mean_hu > 0
        xline(bones{bi}.mean_hu, '-', sprintf('mean=%.0f', bones{bi}.mean_hu), ...
            'Color', colors(bi,:)*0.5, 'LineWidth', 1.5);
    end

    xlim([-400 2000]);
    xlabel('HU');
    ylabel('Voxels');

    if ~isempty(bones{bi}.tag_id)
        tag_str = sprintf(' (tag %d)', bones{bi}.tag_id);
    else
        tag_str = '';
    end
    title(sprintf('Bone %d%s\n%.0f mm^3, mean %.0f HU', ...
        bi, tag_str, bones{bi}.volume_mm3, bones{bi}.mean_hu), 'FontSize', 9);

    % Count any voxels above 2000 HU (possible metal contamination)
    n_metal = sum(hu_vals > 2000);
    if n_metal > 0
        text(0.95, 0.95, sprintf('%d vox > 2000 HU!', n_metal), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'Color', [1 0 0], 'FontWeight', 'bold');
    end
end

sgtitle('HU Distribution per Bone (check for metal contamination)', 'FontSize', 12);

% ========================================================================
%  Figure 4: Mid-slice per-bone close-up
% ========================================================================
if n_bones >= 1
    fig4 = figure('Name', 'Per-Bone Close-ups', 'Color', 'k', ...
        'Position', [200 200 1200 600]);

    n_show = min(n_bones, 5);
    for bi = 1:n_show
        subplot(2, n_show, bi);

        % Find the mid-slice of this bone
        [rr, cc, ss] = ind2sub(size(vol), find(bones{bi}.mask));
        mid_z = round(median(ss));
        bone_slice = bones{bi}.mask(:, :, mid_z);
        ct_slice = vol(:, :, mid_z);

        % Crop to bone region with margin
        r1 = max(1, min(rr)-10); r2 = min(size(vol,1), max(rr)+10);
        c1 = max(1, min(cc)-10); c2 = min(size(vol,2), max(cc)+10);

        ct_crop = ct_slice(r1:r2, c1:c2);
        bone_crop = bone_slice(r1:r2, c1:c2);

        ct_disp = (ct_crop - (-200)) / (1500 - (-200));
        ct_disp = max(0, min(1, ct_disp));
        rgb = repmat(ct_disp, [1 1 3]);

        % Overlay bone
        perim = imdilate(bwperim(bone_crop), strel('disk', 1));
        for ch = 1:3
            layer = rgb(:,:,ch);
            layer(bone_crop) = layer(bone_crop) * 0.5 + colors(bi,ch) * 0.5;
            layer(perim) = colors(bi,ch);
            rgb(:,:,ch) = layer;
        end

        % Overlay markers
        if isfield(sep_result, 'marker_mask')
            mk_crop = sep_result.marker_mask(r1:r2, c1:c2, mid_z);
            if any(mk_crop(:))
                mk_p = imdilate(mk_crop, strel('disk', 1));
                rgb(:,:,1) = max(rgb(:,:,1), double(mk_p));
                rgb(:,:,2) = rgb(:,:,2) .* (1 - double(mk_p)*0.8);
                rgb(:,:,3) = rgb(:,:,3) .* (1 - double(mk_p)*0.8);
            end
        end

        imshow(rgb);
        if ~isempty(bones{bi}.tag_id)
            tag_str = sprintf(' tag%d', bones{bi}.tag_id);
        else
            tag_str = '';
        end
        title(sprintf('Bone %d%s z=%d', bi, tag_str, mid_z), ...
            'Color', 'w', 'FontSize', 9);

        % Bottom row: coronal slice through same bone
        subplot(2, n_show, n_show + bi);
        mid_c = round(median(cc));
        bone_cor = squeeze(bones{bi}.mask(:, mid_c, :));
        ct_cor = squeeze(vol(:, mid_c, :));

        s1 = max(1, min(ss)-10); s2 = min(size(vol,3), max(ss)+10);
        ct_crop_c = ct_cor(r1:r2, s1:s2);
        bone_crop_c = bone_cor(r1:r2, s1:s2);

        ct_disp_c = (ct_crop_c - (-200)) / (1500 - (-200));
        ct_disp_c = max(0, min(1, ct_disp_c));
        rgb_c = repmat(ct_disp_c, [1 1 3]);

        perim_c = imdilate(bwperim(bone_crop_c), strel('disk', 1));
        for ch = 1:3
            layer = rgb_c(:,:,ch);
            layer(bone_crop_c) = layer(bone_crop_c) * 0.5 + colors(bi,ch) * 0.5;
            layer(perim_c) = colors(bi,ch);
            rgb_c(:,:,ch) = layer;
        end

        if isfield(sep_result, 'marker_mask')
            mk_cor = squeeze(sep_result.marker_mask(:, mid_c, :));
            mk_crop_c = mk_cor(r1:r2, s1:s2);
            if any(mk_crop_c(:))
                mk_p_c = imdilate(mk_crop_c, strel('disk', 1));
                rgb_c(:,:,1) = max(rgb_c(:,:,1), double(mk_p_c));
                rgb_c(:,:,2) = rgb_c(:,:,2) .* (1 - double(mk_p_c)*0.8);
                rgb_c(:,:,3) = rgb_c(:,:,3) .* (1 - double(mk_p_c)*0.8);
            end
        end

        imshow(rgb_c);
        title(sprintf('Coronal y=%d', mid_c), 'Color', 'w', 'FontSize', 9);
    end

    sgtitle('Per-Bone Close-ups (top: axial, bottom: coronal)', ...
        'Color', 'w', 'FontSize', 12);
end

% ---- Save figures ----
if isfield(opts, 'OutputDir') && ~isempty(opts.OutputDir)
    outDir = opts.OutputDir;
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    try
        saveas(fig1, fullfile(outDir, 'bone_separation_3d.png'));
        saveas(fig2, fullfile(outDir, 'bone_separation_slices.png'));
        saveas(fig3, fullfile(outDir, 'bone_hu_histograms.png'));
        if exist('fig4', 'var')
            saveas(fig4, fullfile(outDir, 'bone_closeups.png'));
        end
        fprintf('  [Viz] Saved figures to %s\n', outDir);
    catch ME
        warning('Figure save failed: %s', ME.message);
    end
end
end
