function [labels_out, order, means] = relabel_by_hu(labels_in, HU, k, mask)
% RELABEL_BY_HU
% Renumber classes so 1..k correspond to ascending mean HU.
% Returns:
%   labels_out : relabeled uint8 map (0 outside mask)
%   order      : kx1 vector of original labels in ascending HU order
%   means      : kx1 mean HU per ordered class

L = uint16(labels_in);
if nargin < 4 || isempty(mask), mask = (L>0); else, mask = logical(mask); end

means = nan(k,1);
for lbl = 1:k
    msk = (L==lbl) & mask;
    if any(msk(:))
        means(lbl) = mean(HU(msk), 'omitnan');
    end
end
[~, order] = sort(means, 'ascend', 'MissingPlacement','last');

% Build mapping old->new
map = zeros(1, max(k, double(max(L(:)))));
for new = 1:k
    if new <= numel(order) && order(new) >= 1
        map(order(new)) = new;
    end
end

labels_out = zeros(size(L), 'uint8');
in = (L>0) & mask;
labels_out(in) = uint8(map(L(in)));
end
