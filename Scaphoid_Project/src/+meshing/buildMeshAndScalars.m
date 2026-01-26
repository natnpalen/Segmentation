function [mesh, perVertexHU, diag] = buildMeshAndScalars(V, BW, spacing, varargin)
% buildMeshAndScalars  Robust mesh extraction from a binary mask with per-vertex HU sampling.
%
% Usage (backward-compatible):
%   [mesh, perHU] = meshing.buildMeshAndScalars(ds.HU, BW, ds.spacing);
%
% Optional name/value args (all OFF by default to avoid changing behavior):
%   'CloseIfThinMM'       : 0        % If >0, mm-radius closing before iso (helps paper-thin masks)
%   'TrySmoothedFallback' : false    % If true, try isosurface on smooth3(double(BW)) if raw fails
%   'SmoothKernel'        : 'box'    % smooth3 kernel ('box' or 'gaussian')
%   'SmoothSize'          : 3        % smooth3 size (scalar or [r c s])
%   'IsoValue'            : 0.5      % Isovalue for iso extraction
%   'ClampVerticesToGrid' : true     % Clamp mm vertices to valid sampling range
%   'MinVoxels'           : 50       % Fail early if mask smaller than this
%
% Returns:
%   mesh.vertices : [N x 3] in millimeters (LPS-style mm axes: row, col, slice spacing applied)
%   mesh.faces    : [M x 3] triangulation
%   perVertexHU   : [N x 1] trilinear samples from V at mesh vertices
%   diag          : (optional) struct with stats & flags
%
% Notes:
% - This function will NEVER silently invent a surface if the mask is empty.
%   It throws a clear error with diagnostics. Salvage steps are *opt-in*.
% - If you want to allow a gentle "make it watertight" nudge for thin cases,
%   set 'CloseIfThinMM', e.g., to 0.6 (mm).

% ---------------- Options ----------------
opts = struct( ...
  'CloseIfThinMM',        0, ...
  'TrySmoothedFallback',  false, ...
  'SmoothKernel',         'box', ...
  'SmoothSize',           3, ...
  'IsoValue',             0.5, ...
  'ClampVerticesToGrid',  true, ...
  'MinVoxels',            50);
if mod(numel(varargin),2)~=0
  error('Name/value pairs expected.');
end
for k=1:2:numel(varargin)
  name = varargin{k}; val = varargin{k+1};
  if ~isfield(opts, name), error('Unknown option: %s', name); end
  opts.(name) = val;
end

% ------------- Input checks -------------
if ndims(BW)~=3 || ~islogical(BW)
  error('buildMeshAndScalars:badMask', 'BW must be a 3-D logical array.');
end
if ndims(V)~=3 || ~isequal(size(V), size(BW))
  error('buildMeshAndScalars:sizeMismatch', 'V and BW must have identical 3-D size.');
end
if numel(spacing)~=3 || any(~isfinite(spacing)) || any(spacing<=0)
  error('buildMeshAndScalars:badSpacing', 'spacing must be [dr dc ds] in mm, all positive finite.');
end

[R,C,S] = size(BW);
voxCount = nnz(BW);
touchesBorder = false;
if voxCount>0
  border = false(R,C,S); border([1 end],:,:) = true; border(:,[1 end],:) = true; border(:,:,[1 end]) = true;
  touchesBorder = any(BW(border));
end

diag = struct();
diag.size             = [R C S];
diag.voxelCount       = voxCount;
diag.fillFraction     = voxCount / max(1,numel(BW));
diag.touchesBorder    = logical(touchesBorder);
diag.closeAppliedMM   = 0;
diag.usedSmoothFallback = false;
diag.method           = '';
diag.isoValue         = opts.IsoValue;

if voxCount == 0
  error('buildMeshAndScalars:emptyMask', ...
    'Mask is empty (nnz=0). No isosurface can be extracted.');
end
if voxCount < opts.MinVoxels
  error('buildMeshAndScalars:tinyMask', ...
    'Mask too small (nnz=%d < MinVoxels=%d).', voxCount, opts.MinVoxels);
end
if voxCount == numel(BW)
  error('buildMeshAndScalars:fullMask', ...
    'Mask is all-ones. No meaningful surface can be extracted.');
end

% ------------- Optional thin-mask close (mm-aware) -------------
BWwork = BW;
if opts.CloseIfThinMM > 0
  rvox = max(1, round(opts.CloseIfThinMM / max(spacing(:))));
  SE = strel('sphere', rvox);
  BWwork = imclose(BWwork, SE);
  diag.closeAppliedMM = opts.CloseIfThinMM;
  if ~any(BWwork(:))
    error('buildMeshAndScalars:closeWipedMask', ...
      'CloseIfThinMM removed the mask unexpectedly. Disable or reduce radius.');
  end
end

% ------------- Primary ISO extraction -------------
[faces, vertices] = try_extract_iso(BWwork, opts.IsoValue);
diag.method = 'raw';

% ------------- Fallback: smoothed mask -------------
if (isempty(vertices) || isempty(faces)) && opts.TrySmoothedFallback
  BWdbl = double(BWwork);
  sz = opts.SmoothSize;
  if isscalar(sz), sz = [sz sz sz]; end
  try
    Vs = smooth3(BWdbl, opts.SmoothKernel, sz);
  catch
    % Some MATLABs require [r c s] ints; enforce
    sz = max(1, round(sz));
    Vs = smooth3(BWdbl, opts.SmoothKernel, sz);
  end
  [faces, vertices] = try_extract_iso(Vs, opts.IsoValue);
  diag.method = 'smooth-fallback';
  diag.usedSmoothFallback = true;
end

% If still empty: fail with diagnostics
if isempty(vertices) || isempty(faces)
  error('buildMeshAndScalars:isosurfaceFailed', ...
    ['Failed to extract an isosurface from the mask. Diagnostics:\n' ...
     '  nnz(BW)=%d (%.6f fill)\n' ...
     '  touchesBorder=%d\n' ...
     '  CloseIfThinMM=%.3f (applied=%.3f)\n' ...
     '  TrySmoothedFallback=%d\n'], ...
     voxCount, diag.fillFraction, diag.touchesBorder, opts.CloseIfThinMM, diag.closeAppliedMM, opts.TrySmoothedFallback);
end

% ------------- Convert iso coords to mm (axis reorder) -------------
% isosurface returns vertices in voxel index space with coords:
%   x = columns, y = rows, z = slices, each 1-based-ish
verts_vox = vertices - 1; % make it 0-based
verts_mm  = [verts_vox(:,2)*spacing(1), ...
             verts_vox(:,1)*spacing(2), ...
             verts_vox(:,3)*spacing(3)];

% ------------- Optional clamp to valid mm domain -------------
if opts.ClampVerticesToGrid
  limR = (R-1)*spacing(1);
  limC = (C-1)*spacing(2);
  limS = (S-1)*spacing(3);
  verts_mm(:,1) = max(0, min(limR, verts_mm(:,1)));
  verts_mm(:,2) = max(0, min(limC, verts_mm(:,2)));
  verts_mm(:,3) = max(0, min(limS, verts_mm(:,3)));
end

% ------------- Per-vertex HU (trilinear) -------------
[Xmm,Ymm,Zmm] = ndgrid( (0:R-1)*spacing(1), (0:C-1)*spacing(2), (0:S-1)*spacing(3) );
perVertexHU = interp3( Ymm, Xmm, Zmm, double(V), ...
                       verts_mm(:,2), verts_mm(:,1), verts_mm(:,3), ...
                       'linear', NaN);

mesh = struct('vertices', verts_mm, 'faces', faces);

% ------------- Local helpers -------------
function [F, Vtx] = try_extract_iso(vol, iso)
  F = []; Vtx = [];
  % Prefer Medical Imaging Toolbox if present
  if exist('extractIsosurface','file')==2
    try
      fv = extractIsosurface(vol, iso);
      if isa(fv, 'triangulation')
        F   = fv.ConnectivityList;
        Vtx = fv.Points;
        if ~isempty(F) && ~isempty(Vtx), return; end
      elseif isstruct(fv) && isfield(fv,'faces') && isfield(fv,'vertices')
        F   = fv.faces;
        Vtx = fv.vertices;
        if ~isempty(F) && ~isempty(Vtx), return; end
      else
        % Some releases support [F,V] outputs
        try
          [F,Vtx] = extractIsosurface(vol, iso);
          if ~isempty(F) && ~isempty(Vtx), return; end
        catch, end
      end
    catch
      % fall through to base isosurface
    end
  end
  % Base MATLAB isosurface
  try
    fv2 = isosurface(vol, iso);      % struct output
    F   = fv2.faces;
    Vtx = fv2.vertices;
  catch
    try
      [F,Vtx] = isosurface(vol, iso); % 2-output legacy form
    catch
      F = []; Vtx = [];
    end
  end
end
end
