function show_mesh_3d(mesh, hu)
% Simple 3-D viewer for the output mesh; colors by HU if provided
figure('Name','Scaphoid Mesh','Color','w');
p = patch('Faces',mesh.faces, 'Vertices',mesh.vertices, ...
       'FaceColor',[0.85 0.85 0.9], 'EdgeColor','none', 'FaceAlpha',1);
if nargin>1 && ~isempty(hu) && numel(hu)==size(mesh.vertices,1)
 set(p,'FaceVertexCData', hu, 'FaceColor','interp');
 colorbar; ylabel(colorbar,'HU');
end
axis equal off vis3d
camorbit(30,10); camlight headlight; lighting gouraud
end
