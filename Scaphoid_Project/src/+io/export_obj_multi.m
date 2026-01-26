function export_obj_multi(meshes, names, outObjFile, varargin)
% EXPORT_OBJ_MULTI  Write multiple meshes into one OBJ with per-object names & MTL.
% meshes    : cell of {V,F}
% names     : cell of strings for 'o <name>' blocks
% outObjFile: path/to/file.obj
% Options:
%   'MaterialNames' : cellstr same length as meshes (default: derived from names)
%   'DiffuseColors' : Nx3 double [0..1] per mesh (optional; for MTL)

p = inputParser;
addParameter(p,'MaterialNames',{},@(x)iscell(x));
addParameter(p,'DiffuseColors',[],@(x)isnumeric(x) && (isempty(x) || size(x,2)==3));
parse(p,varargin{:});
matNames = p.Results.MaterialNames;
Kd = p.Results.DiffuseColors;

n = numel(meshes);
if isempty(matNames), matNames = names; end
if isempty(Kd), Kd = repmat([0.8 0.8 0.8], n,1); end

[objDir, base, ~] = fileparts(outObjFile);
if objDir==""; objDir = pwd; end
outMtlFile = fullfile(objDir, base + ".mtl");
mtlRel = string(base + ".mtl");

% --- Write MTL ---
fid = fopen(outMtlFile, 'w');
assert(fid>0, 'Cannot open %s', outMtlFile);
fprintf(fid, '# Multi-material for %s\n', base);
for i=1:n
    name = sanitize(matNames{i});
    fprintf(fid, 'newmtl %s\n', name);
    fprintf(fid, 'Kd %.6f %.6f %.6f\n', Kd(i,1), Kd(i,2), Kd(i,3));
    fprintf(fid, 'Ka 0 0 0\nKs 0 0 0\nNs 1\nillum 1\n\n');
end
fclose(fid);

% --- Write OBJ ---
fo = fopen(outObjFile, 'w');
assert(fo>0, 'Cannot open %s', outObjFile);
fprintf(fo, '# Multi-object OBJ\nmtllib %s\n', mtlRel);

vOffset = 0;
for i=1:n
    m = meshes{i};
    if isempty(m) || ~isfield(m,'vertices') || isempty(m.vertices) || ~isfield(m,'faces') || isempty(m.faces)
        continue;
    end
    V = double(m.vertices); F = double(m.faces);
    name = sanitize(names{i});
    mat  = sanitize(matNames{i});

    fprintf(fo, '\no %s\nusemtl %s\n', name, mat);

    % Vertices
    for vi=1:size(V,1)
        fprintf(fo, 'v %.6f %.6f %.6f\n', V(vi,1), V(vi,2), V(vi,3));
    end
    % Faces (1-based in OBJ; add offset)
    for fi=1:size(F,1)
        fprintf(fo, 'f %d %d %d\n', F(fi,1)+vOffset, F(fi,2)+vOffset, F(fi,3)+vOffset);
    end
    vOffset = vOffset + size(V,1);
end

fclose(fo);
fprintf('OBJ+MTL written: %s, %s\n', outObjFile, outMtlFile);
end

function s = sanitize(s)
s = regexprep(string(s), '[^\w\-\.\(\)]+','_');
end
