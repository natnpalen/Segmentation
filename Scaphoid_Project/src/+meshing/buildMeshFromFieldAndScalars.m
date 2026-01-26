function [mesh, perVertexHU] = buildMeshFromFieldAndScalars(field, isoval, spacing, HUvolume)
% Build mesh from a scalar field (double), sampling HU at vertices.
field    = double(field);
HUvolume = double(HUvolume);

% Quick sanity: ensure level exists
fmin = min(field(:)); fmax = max(field(:));
if ~(fmin <= isoval && isoval <= fmax)
    error('buildMeshFromFieldAndScalars:isosurfaceFailed', ...
          'Field has no %g-level crossing (min=%g, max=%g).', isoval, fmin, fmax);
end

fv = isosurface(field, isoval);     % base MATLAB path for consistency
faces = fv.faces; vertices = fv.vertices;

% Convert iso coords [x=cols,y=rows,z=slices] → mm with axis reorder
verts_vox = vertices - 1;
verts_mm  = [verts_vox(:,2)*spacing(1), ...
             verts_vox(:,1)*spacing(2), ...
             verts_vox(:,3)*spacing(3)];

% Per-vertex HU via trilinear interpolation in mm coords
[R,C,S] = size(HUvolume);
[Xmm,Ymm,Zmm] = ndgrid( (0:R-1)*spacing(1), (0:C-1)*spacing(2), (0:S-1)*spacing(3) );
perVertexHU = interp3( Ymm, Xmm, Zmm, HUvolume, ...
                       verts_mm(:,2), verts_mm(:,1), verts_mm(:,3), ...
                       'linear', NaN);

mesh = struct('vertices', verts_mm, 'faces', faces);
end
