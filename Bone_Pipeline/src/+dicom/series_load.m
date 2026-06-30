function ds = series_load(folder, varargin)
% SERIES_LOAD  Load a DICOM CT series into a structured dataset.
%
%   ds = dicom.series_load(folder)
%   ds = dicom.series_load(folder, 'TargetIsoMM', 0.5, 'Smoothing', true)
%
% Reads all DICOM files in the folder, sorts by slice position, applies
% rescale slope/intercept to produce HU values, and returns a struct with
% the volume, spacing, orientation, and coordinate transforms.
%
% Uses direct file enumeration (robust for scanner-exported DICOM without
% standard extensions or complete metadata in every file).

o = struct('TargetIsoMM', [], 'Smoothing', false);
o = utils.parse_opts(o, varargin{:});

% ---- Find valid DICOM files ----
files = list_dicom_files(folder);
fprintf('[Load] Found %d DICOM files in %s\n', numel(files), folder);
if isempty(files)
    error('No readable DICOM files found under: %s', folder);
end

% ---- Read all slices ----
[Vraw, infosCell] = read_stack(files);
fprintf('[Load] Volume loaded: size %s\n', mat2str(size(Vraw)));

% ---- Normalize shape: [R C 1 S] -> [R C S] ----
Vraw = squeeze_volume(Vraw);

% ---- Convert to HU ----
HU = apply_rescale(double(Vraw), infosCell);

% ---- Orientation, ordering, spacing ----
[dr, dc, ds_sp, origin, dir_row, dir_col, dir_slice, HU] = ...
    orient_and_order(HU, infosCell);

% ---- Optional smoothing ----
if o.Smoothing
    G = imgradient3(HU);
    edgeMask = G > prctile(G(:), 95);
    HU(~edgeMask) = imdiffusefilt(HU(~edgeMask), 'ConductionMethod', 'quadratic', ...
        'NumberOfIterations', 2, 'GradientThreshold', 25);
end

% ---- Optional isotropic resampling ----
if ~isempty(o.TargetIsoMM)
    HU = squeeze_volume(HU);
    [HU, dr, dc, ds_sp] = resample_isotropic(HU, dr, dc, ds_sp, o.TargetIsoMM);
end

% ---- Build output ----
[R, C, S] = size(HU);

voxelToWorld = @(i_row, i_col, i_slice) voxel_to_world( ...
    i_row, i_col, i_slice, origin, dir_row, dir_col, dir_slice, dr, dc, ds_sp);
worldToVoxel = @(xyz) world_to_voxel( ...
    xyz, origin, dir_row, dir_col, dir_slice, dr, dc, ds_sp);

ds = struct();
ds.HU          = HU;
ds.spacing     = [dr dc ds_sp];
ds.origin      = origin(:).';
ds.dir_row     = dir_row(:).';
ds.dir_col     = dir_col(:).';
ds.dir_slice   = dir_slice(:).';
ds.size        = [R C S];
ds.infos       = infosCell;
ds.voxelToWorld = voxelToWorld;
ds.worldToVoxel = worldToVoxel;
ds.M_voxToLPS  = [dir_col*dc, dir_row*dr, dir_slice*ds_sp, origin(:); 0 0 0 1];
fprintf('[Load] OK: %dx%dx%d, spacing [%.3f %.3f %.3f] mm\n', R, C, S, dr, dc, ds_sp);
end


% =========================================================================
%  FILE DISCOVERY
% =========================================================================
function filesList = list_dicom_files(folder)
    if ~isfolder(folder), error('Folder not found: %s', folder); end
    listing = dir(fullfile(folder, '**', '*'));
    files = listing(~[listing.isdir]);
    paths = fullfile({files.folder}, {files.name});
    keep = false(size(paths));
    for i = 1:numel(paths)
        try
            d = dir(paths{i});
            if d.bytes < 512, continue; end
            dicominfo(paths{i});
            keep(i) = true;
        catch
        end
    end
    filesList = paths(keep);
end


% =========================================================================
%  MANUAL STACK READ
% =========================================================================
function [V, infosCell] = read_stack(filesList)
    meta = struct('fname', [], 'InstanceNumber', [], 'IPP', [], 'IOP', []);
    M = repmat(meta, 0, 1);
    for i = 1:numel(filesList)
        try
            info = dicominfo(filesList{i});
            rec.fname = filesList{i};
            rec.InstanceNumber = get_field(info, 'InstanceNumber', i);
            rec.IOP = get_field(info, 'ImageOrientationPatient', [1 0 0 0 1 0]);
            if isfield(info, 'ImagePositionPatient')
                rec.IPP = double(info.ImagePositionPatient(:)).';
            else
                rec.IPP = [NaN NaN NaN];
            end
            M(end+1, 1) = rec; %#ok<AGROW>
        catch
        end
    end
    if isempty(M), error('No readable DICOM headers.'); end

    % Sort by slice position along the normal
    IOP0 = M(1).IOP;
    drow = IOP0(1:3); drow = drow(:) / max(norm(drow), eps);
    dcol = IOP0(4:6); dcol = dcol(:) / max(norm(dcol), eps);
    nrm = cross(drow, dcol); nrm = nrm / max(norm(nrm), eps);

    hasIPP = all(isfinite(M(1).IPP));
    if hasIPP
        projs = cellfun(@(p) dot(p, nrm), {M.IPP});
        [~, ord] = sort(projs, 'ascend');
    else
        [~, ord] = sort([M.InstanceNumber], 'ascend');
    end
    M = M(ord);

    info0 = dicominfo(M(1).fname);
    R = double(get_field(info0, 'Rows', []));
    C = double(get_field(info0, 'Columns', []));
    S = numel(M);

    infosCell = cell(S, 1);
    V = zeros(R, C, S, 'double');
    for k = 1:S
        infosCell{k} = dicominfo(M(k).fname);
        V(:,:,k) = double(dicomread(infosCell{k}));
    end
end


% =========================================================================
%  RESCALE SLOPE/INTERCEPT
% =========================================================================
function HU = apply_rescale(V, infosCell)
    S = size(V, 3);
    if isempty(infosCell)
        HU = V;
        return;
    end
    if numel(infosCell) == 1
        [s, b] = slope_intercept(infosCell{1});
        HU = s .* V + b;
        return;
    end

    slopes = ones(S, 1);
    intercepts = zeros(S, 1);
    n = min(S, numel(infosCell));
    for k = 1:n
        [slopes(k), intercepts(k)] = slope_intercept(infosCell{k});
    end
    if n < S
        slopes(n+1:S) = slopes(n);
        intercepts(n+1:S) = intercepts(n);
    end

    HU = zeros(size(V));
    for k = 1:S
        HU(:,:,k) = slopes(k) .* V(:,:,k) + intercepts(k);
    end
end


function [s, b] = slope_intercept(info)
    s = 1; b = 0;
    if isempty(info), return; end
    if isfield(info, 'RescaleSlope'), s = double(info.RescaleSlope); end
    if isfield(info, 'RescaleIntercept'), b = double(info.RescaleIntercept); end
    if ~isfinite(s), s = 1; end
    if ~isfinite(b), b = 0; end
end


% =========================================================================
%  ORIENTATION AND ORDERING
% =========================================================================
function [dr, dc, ds, origin, dir_row, dir_col, dir_slice, HU] = orient_and_order(HU, infosCell)
    dr = 1; dc = 1; ds = 1;
    origin = [0; 0; 0];
    dir_row = [1; 0; 0]; dir_col = [0; 1; 0]; dir_slice = [0; 0; 1];

    Svol = size(HU, 3);
    if isempty(infosCell), return; end

    info0 = infosCell{1};
    if isfield(info0, 'PixelSpacing') && numel(info0.PixelSpacing) >= 2
        dr = double(info0.PixelSpacing(1));
        dc = double(info0.PixelSpacing(2));
    end
    if isfield(info0, 'ImageOrientationPatient') && numel(info0.ImageOrientationPatient) >= 6
        IOP = double(info0.ImageOrientationPatient(:));
        dir_row = IOP(1:3); dir_row = dir_row / max(norm(dir_row), eps);
        dir_col = IOP(4:6); dir_col = dir_col / max(norm(dir_col), eps);
        dir_slice = cross(dir_row, dir_col);
        dir_slice = dir_slice / max(norm(dir_slice), eps);
    end

    if isfield(info0, 'SpacingBetweenSlices') && isfinite(info0.SpacingBetweenSlices)
        ds = double(info0.SpacingBetweenSlices);
    elseif isfield(info0, 'SliceThickness') && isfinite(info0.SliceThickness)
        ds = double(info0.SliceThickness);
    end

    Smeta = numel(infosCell);
    if Smeta == Svol
        [IPP_all, maskIPP] = collect_ipp(infosCell);
        nIPP = sum(maskIPP);
        if nIPP >= 2
            proj = IPP_all(maskIPP,:) * dir_slice;
            dproj = abs(diff(sort(proj, 'ascend')));
            if ~isempty(dproj), ds = median(dproj(~isnan(dproj))); end
        end
        if nIPP == Svol
            [~, ord] = sort(IPP_all * dir_slice, 'ascend');
            HU = HU(:,:,ord);
            infosCell = infosCell(ord); %#ok<NASGU>
            origin = IPP_all(ord(1),:).';
            return;
        end
        if nIPP >= 1
            firstIdx = find(maskIPP, 1, 'first');
            origin = IPP_all(firstIdx,:).';
        end
    else
        if isfield(info0, 'ImagePositionPatient') && numel(info0.ImagePositionPatient) == 3
            origin = double(info0.ImagePositionPatient(:));
        end
    end
end


function [IPP_all, mask] = collect_ipp(infosCell)
    S = numel(infosCell);
    IPP_all = NaN(S, 3); mask = false(S, 1);
    for k = 1:S
        if isfield(infosCell{k}, 'ImagePositionPatient')
            v = double(infosCell{k}.ImagePositionPatient(:)).';
            if numel(v) == 3 && all(isfinite(v))
                IPP_all(k,:) = v; mask(k) = true;
            end
        end
    end
end


% =========================================================================
%  ISOTROPIC RESAMPLING
% =========================================================================
function [HU_iso, dr, dc, ds] = resample_isotropic(HU, dr, dc, ds, iso_mm)
    HU = squeeze(HU);
    [R, C, S] = size(HU);
    F = griddedInterpolant({1:R, 1:C, 1:S}, double(HU), 'linear', 'nearest');

    newR = max(1, round((R-1)*dr/iso_mm) + 1);
    newC = max(1, round((C-1)*dc/iso_mm) + 1);
    newS = max(1, round((S-1)*ds/iso_mm) + 1);

    rq = ((0:newR-1)*iso_mm)/dr + 1;
    cq = ((0:newC-1)*iso_mm)/dc + 1;
    sq = ((0:newS-1)*iso_mm)/ds + 1;
    [RR, CC, SS] = ndgrid(rq, cq, sq);
    HU_iso = F(RR, CC, SS);
    dr = iso_mm; dc = iso_mm; ds = iso_mm;
end


% =========================================================================
%  COORDINATE TRANSFORMS
% =========================================================================
function [x, y, z] = voxel_to_world(ir, ic, is, origin, dr_v, dc_v, ds_v, dr, dc, dsl)
    ir = double(ir); ic = double(ic); is = double(is);
    d_r = (ir-1).*dr; d_c = (ic-1).*dc; d_s = (is-1).*dsl;
    dx = d_c.*dc_v(1) + d_r.*dr_v(1) + d_s.*ds_v(1);
    dy = d_c.*dc_v(2) + d_r.*dr_v(2) + d_s.*ds_v(2);
    dz = d_c.*dc_v(3) + d_r.*dr_v(3) + d_s.*ds_v(3);
    x = origin(1)+dx; y = origin(2)+dy; z = origin(3)+dz;
end


function ijk = world_to_voxel(xyz, origin, dr_v, dc_v, ds_v, dr, dc, dsl)
    d = bsxfun(@minus, xyz, origin(:).');
    ijk = [d*dr_v/dr + 1, d*dc_v/dc + 1, d*ds_v/dsl + 1];
end


% =========================================================================
%  UTILITIES
% =========================================================================
function V = squeeze_volume(V)
    if ndims(V) == 4
        if size(V,3) == 1, V = squeeze(V); end
        if size(V,4) == 1, V = squeeze(V); end
    end
end


function val = get_field(s, field, def)
    if isfield(s, field), val = s.(field); else, val = def; end
end
