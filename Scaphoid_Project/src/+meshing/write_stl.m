function write_stl(filename, mesh)
% writeSTL  Write a triangulated surface to ASCII STL.
% mesh must have fields .vertices (Nx3) and .faces (Mx3)
V = mesh.vertices; F = mesh.faces;
fid = fopen(filename,'w');
if fid < 0, error('Cannot open %s for writing.', filename); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'solid scaphoid\n');
for i=1:size(F,1)
 tri = F(i,:);
 v1 = V(tri(1),:); v2 = V(tri(2),:); v3 = V(tri(3),:);
 n = cross(v2 - v1, v3 - v1);
 nn = norm(n); if nn>0, n = n/nn; else, n = [0 0 0]; end
 fprintf(fid, '  facet normal %.6g %.6g %.6g\n', n);
 fprintf(fid, '    outer loop\n');
 fprintf(fid, '      vertex %.6g %.6g %.6g\n', v1);
 fprintf(fid, '      vertex %.6g %.6g %.6g\n', v2);
 fprintf(fid, '      vertex %.6g %.6g %.6g\n', v3);
 fprintf(fid, '    endloop\n');
 fprintf(fid, '  endfacet\n');
end
fprintf(fid, 'endsolid scaphoid\n');
end
