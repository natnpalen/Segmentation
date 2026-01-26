function hires_data = upsample_data(lowres_data, upsamplingFactor)
% UPSAMPLE_DATA  Upsample HU volume and mask using imresize3 (Image Processing Toolbox).
%
%   Robust masked upsampling:
%     - Resizes HU with a normalized (num/den) scheme to avoid NaN/edge bleed.
%     - Resizes mask with 'nearest' to preserve crisp boundaries.
%     - Sets HU to NaN outside the upsampled mask.
%
%   SYNTAX:
%       hires_data = partition.upsample_data(lowres_data, upsamplingFactor)
%
%   INPUTS:
%       lowres_data (struct) with fields:
%           .HU            (3D numeric): Hounsfield Units. Outside mask may be NaN or -3000.
%           .scaphoid_mask (3D logical): object mask (same size as HU).
%           .ds            (struct)     : must have ds.spacing = [dr dc ds] (mm/voxel).
%       upsamplingFactor (scalar > 0): integer factor to increase resolution.
%
%   OUTPUT:
%       hires_data (struct):
%           .HU   (3D double)  : upsampled HU, NaN outside mask.
%           .mask (3D logical) : upsampled mask.
%           .ds   (struct)     : updated dataset with spacing/size (+ geometry passthrough).
%
%   Notes:
%     - Uses 'cubic' for HU and mask-weight (smooth), 'nearest' for final mask.
%     - Avoids spline overshoot and NaN propagation by masked normalization.

% ---- Unpack & validate ---------------------------------------------------
HU_lowres   = double(lowres_data.HU);
mask_lowres = logical(lowres_data.scaphoid_mask);
ds_lowres   = lowres_data.ds;

if ~(isnumeric(upsamplingFactor) && isscalar(upsamplingFactor) && upsamplingFactor > 0)
    error('upsamplingFactor must be a positive scalar.');
end
if ~isfield(ds_lowres, 'spacing') || numel(ds_lowres.spacing) ~= 3
    error('lowres_data.ds.spacing must be a 3-element vector [dr dc ds] in mm/voxel.');
end

[R, C, S] = size(HU_lowres);
spacing    = double(ds_lowres.spacing(:)).';   % [dr dc ds] mm/voxel

% Normalize background:
%   - If upstream stored -3000 outside, convert that to NaN for clarity here.
%   - We'll use masked normalization, so we don't rely on NaNs during resize.
HU_lowres(HU_lowres <= -3000 & ~mask_lowres) = NaN;

% ---- Compute hi-res size -------------------------------------------------
f = double(upsamplingFactor);
R_hires = (R-1) * f + 1;
C_hires = (C-1) * f + 1;
S_hires = (S-1) * f + 1;
targetSz = [R_hires, C_hires, S_hires];

% --- NaN-safe masked resize (important!) ---
% Ensure no NaNs go into imresize3:
HU_filled = HU_lowres;
HU_filled(~mask_lowres | isnan(HU_filled)) = 0;    % ← zero out outside & any NaNs

mask_weight_lowres = double(mask_lowres);          % 0/1 weights

% Smooth numerator/denominator independently, *without* NaNs:
num = imresize3(HU_filled .* mask_weight_lowres, targetSz, 'cubic');
den = imresize3(mask_weight_lowres,              targetSz, 'cubic');

% Normalize; where den≈0 (outside), result will be 0 for now—we'll set to NaN below.
HU_hires = num ./ max(den, eps);

% ---- Upsample mask (final binary) with 'nearest' ------------------------
mask_hires = imresize3(mask_lowres, targetSz, 'nearest') > 0.5;

% Force HU outside mask to NaN
HU_hires(~mask_hires) = NaN;

% ---- Build ds_hires (spacing & geometry passthrough) --------------------
ds_hires = struct();
ds_hires.spacing = spacing / f;                % mm/voxel at hi-res
ds_hires.size    = targetSz;

% Pass through helpful geometry fields if present
pass = {'origin','dir_row','dir_col','dir_slice','M_voxToLPS','PatientName','SeriesDescription'};
for k = 1:numel(pass)
    fld = pass{k};
    if isfield(ds_lowres, fld)
        ds_hires.(fld) = ds_lowres.(fld);
    end
end
% Clarify description if present
if isfield(ds_hires,'SeriesDescription')
    ds_hires.SeriesDescription = sprintf('%s (Upsampled %dx)', string(ds_hires.SeriesDescription), f);
else
    ds_hires.SeriesDescription = sprintf('Upsampled Scaphoid (%dx)', f);
end

% ---- Package outputs -----------------------------------------------------
hires_data = struct();
hires_data.HU   = HU_hires;
hires_data.mask = mask_hires;
hires_data.ds   = ds_hires;

% ---- Final sanity checks -------------------------------------------------
assert(all(isfinite(hires_data.ds.spacing)) && all(hires_data.ds.spacing > 0), 'Invalid spacing in ds_hires.');
assert(isequal(size(hires_data.HU), size(hires_data.mask)), 'HU and mask sizes differ after upsampling.');

end
