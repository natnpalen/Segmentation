function [HU, mask, ds] = load_masked_hu(outputDir, fallbackMask, fallbackDS)
% LOAD_MASKED_HU  Robustly load HU and mask aligned in the same grid.
% Prefers saved NIfTI created by your pipeline; otherwise falls back to in-memory.

nii_path = fullfile(outputDir, 'scaphoid_masked_hu.nii.gz');
hu_mat   = fullfile(outputDir, 'HU_masked.mat');
maskinfo = fullfile(outputDir, 'mask_info.mat');

if exist(nii_path, 'file')
    HU = double(niftiread(nii_path));
    HU(HU <= -3000) = NaN;            % restore NaNs that were baked for outside
    % Mask comes from fallback (authoritative segmentation)
    mask = logical(fallbackMask);
    ds   = fallbackDS;
elseif exist(hu_mat, 'file') && exist(maskinfo, 'file')
    S = load(hu_mat);      % HU_masked
    MI = load(maskinfo);   % mask_info
    HU = double(S.HU_masked);
    mask = logical(fallbackMask);
    ds   = fallbackDS;
else
    % Fallback: use ds.HU + mask from pipeline output
    ds   = fallbackDS;
    mask = logical(fallbackMask);
    HU   = double(ds.HU);
    HU(~mask) = NaN;
end
end
