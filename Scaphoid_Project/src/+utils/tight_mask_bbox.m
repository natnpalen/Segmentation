function [bbox, Msub] = tight_mask_bbox(M, marginVox)
% tightMaskBBox: return tight index ranges (with margin) around a 3-D logical mask.
% bbox has fields .r, .c, .s (each a vector of indices).
if nargin<2, marginVox = 0; end
[R,C,S] = size(M);
[idxR, idxC, idxS] = ind2sub([R,C,S], find(M));
if isempty(idxR)
  % empty mask: return full volume to avoid errors
  bbox = struct('r', 1:R, 'c', 1:C, 's', 1:S);
  Msub = M;
  return;
end
r1 = max(1, min(idxR) - marginVox); r2 = min(R, max(idxR) + marginVox);
c1 = max(1, min(idxC) - marginVox); c2 = min(C, max(idxC) + marginVox);
s1 = max(1, min(idxS) - marginVox); s2 = min(S, max(idxS) + marginVox);
bbox = struct('r', r1:r2, 'c', c1:c2, 's', s1:s2);
Msub = M(bbox.r, bbox.c, bbox.s);
end