function labels_out = remove_islands_and_majority(labels_in, mask, k, island_min_mm3, majority_iters, nbrhood_size_vox, spacing)
% REMOVE_ISLANDS_AND_MAJORITY
% 1) Remove small 26-connected components per label (size in mm^3)
% 2) Majority vote to refill 0-labeled voxels INSIDE the mask
% 3) Ensure 0 outside mask

Z = uint8(labels_in); Z(~mask)=0;

voxel_mm3 = prod(double(spacing(:)'));
min_vox   = max(1, round(double(island_min_mm3) / voxel_mm3));

% Remove tiny components per label → set to 0 for reassignment
for lbl = 1:k
    bw = (Z == lbl);
    if ~any(bw(:)), continue; end
    CC = bwconncomp(bw, 26);
    if CC.NumObjects == 0, continue; end
    comp_sizes = cellfun(@numel, CC.PixelIdxList);
    small_ids  = find(comp_sizes < min_vox);
    for s = 1:numel(small_ids)
        Z(CC.PixelIdxList{small_ids(s)}) = 0;
    end
end

% Majority vote passes
ker = ones(nbrhood_size_vox);
for it = 1:majority_iters
    to_fill = (Z == 0) & mask;
    if ~any(to_fill(:)), break; end
    counts = zeros([size(Z) k], 'double');
    for lbl = 1:k
        counts(:,:,:,lbl) = convn(double(Z == lbl), ker, 'same');
    end
    [~, winner] = max(counts, [], 4);
    Z(to_fill) = uint8(winner(to_fill));
end

Z(~mask) = 0;
labels_out = uint8(Z);
end
