function labels_s = smooth_labels_masked(labels_in, mask, k, neigh_mm, spacing)
% SMOOTH_LABELS_MASKED  Mask-aware channel-wise median smoothing on labels.
% Ensures: 0 outside mask, restores raw label on any in-mask holes created.

labels_in  = uint8(labels_in);
labels_in(~mask) = 0;

win = density.compute_spacing_window(neigh_mm, spacing);

% Channel-wise 3D median → argmax
best = zeros(size(labels_in), 'uint8');
acc  = zeros(size(labels_in), 'uint8');
for lab = 1:k
    ch   = uint8(labels_in == lab);
    ch_s = medfilt3(ch, win);
    mask_better = ch_s > acc;
    best(mask_better) = uint8(lab);
    acc(mask_better)  = ch_s(mask_better);
end

labels_s = zeros(size(labels_in), 'uint8');
labels_s(mask) = best(mask);

% Restore raw where smoothing produced zeros inside mask
holes = mask & (labels_s == 0);
if any(holes(:))
    labels_s(holes) = labels_in(holes);
end
labels_s(~mask) = 0;
end
