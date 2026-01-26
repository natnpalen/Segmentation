function BW1 = keep_largest_component_3d(BW)
% keepLargestComponent3D  Keep only the largest 26-connected component.
CC = bwconncomp(logical(BW), 26);
if CC.NumObjects <= 1
  BW1 = logical(BW);
  return;
end
[~, iMax] = max(cellfun(@numel, CC.PixelIdxList));
BW1 = false(size(BW));
BW1(CC.PixelIdxList{iMax}) = true;
end