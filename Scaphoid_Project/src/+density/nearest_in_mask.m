function vals = nearest_in_mask(HU, mask, q)
% nearest value inside mask (q: query logical)
D = bwdist(mask);                     %#ok<NASGU> % (not used, but left if you want speeds up with knnsearch)
[idx_r, idx_c, idx_s] = ind2sub(size(mask), find(mask));
P = [idx_r, idx_c, idx_s];
[qr, qc, qs] = ind2sub(size(q), find(q));
Q = [qr, qc, qs];
% Brute nearest neighbor (small rim only): 
Mdl = createns(P);                    % Statistics Toolbox
nn = knnsearch(Mdl, Q, 'K', 1);
vals = HU(sub2ind(size(HU), P(nn,1), P(nn,2), P(nn,3)));
end
