function fv = isosurface_mm(field, isoval, spacing)
% Base-MATLAB isosurface, then convert voxel coords → mm.
fv0 = isosurface(field, isoval);           % uses indices (x=cols,y=rows,z=slices)
verts_vox = fv0.vertices - 1;
verts_mm  = [verts_vox(:,2)*spacing(1), ...
             verts_vox(:,1)*spacing(2), ...
             verts_vox(:,3)*spacing(3)];
fv = struct('faces', fv0.faces, 'vertices', verts_mm);
end
