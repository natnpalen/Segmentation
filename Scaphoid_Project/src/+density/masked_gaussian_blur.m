function Vb = masked_gaussian_blur(V, mask, sigma_mm, spacing)
% MASKED_GAUSSIAN_BLUR  Blur V inside mask in a NaN-safe way; keep outside as NaN.
% V: numeric 3D volume (NaNs outside okay), mask: logical, sigma_mm: scalar or [x y z].

if ~any(sigma_mm), Vb = V; Vb(~mask) = NaN; return; end
sigma_mm = double(sigma_mm);
if isscalar(sigma_mm), sigma_mm = [sigma_mm sigma_mm sigma_mm]; end
sig_vox  = max(eps, sigma_mm ./ double(spacing(:)'));

Vd   = double(V); Vd(~mask) = 0;
Mb   = imgaussfilt3(double(mask), sig_vox);
Vb_n = imgaussfilt3(Vd,            sig_vox);
Vb   = Vb_n ./ max(Mb, eps);
Vb(~mask) = NaN;
end
