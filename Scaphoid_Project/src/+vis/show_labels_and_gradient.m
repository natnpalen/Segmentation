function show_labels_and_gradient(labels_kmeans, gradient_map, ds, k, mask)
% SHOW_LABELS_AND_GRADIENT  One-call viewer for labels & binned gradient.

cmap_kmeans = [0.85 0.85 0.85; parula(k)];
alpha_mask  = double(mask);

scroll_viewer(labels_kmeans, ds, 'K-Means Zones (Smoothed)', ...
    cmap_kmeans, 'Zone ID', 'clim', [0 k], 'alpha', alpha_mask);

cmap_grad = [0 0 0; hot(256)];
grad_disp = zeros(size(gradient_map), 'uint16');
inb = mask & isfinite(gradient_map);
grad_disp(inb) = 1 + floor(255 * gradient_map(inb));

scroll_viewer(grad_disp, ds, 'Continuous Density Gradient (Binned)', ...
    cmap_grad, 'Relative Density (bins)', 'clim', [0 256]);
end
