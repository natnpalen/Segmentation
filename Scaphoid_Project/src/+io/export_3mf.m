function export_3mf(meshes, names, outFile)
% EXPORT_3MF  Stream multiple meshes to a 3MF (UTF-8), Cura/Bambu‑friendly.
% FIXED: preserves triangle winding; de‑dupes faces without altering order;
%        safer quad triangulation; robust cleaning; assembly wrapper.
%
% Inputs
%   meshes : cell array of structs with fields .vertices (Nx3), .faces (MxK)
%   names  : cell array of part names (same length as meshes)
%   outFile: string path to .3mf
%
% Notes
% - This version deliberately DOES NOT sort vertex indices inside faces.
%   Slicers rely on consistent winding; changing order breaks orientation
%   and can lead to tiny/negative volumes and huge non‑manifold counts.
% - Duplicate faces are removed using an order‑independent key, while
%   keeping the ORIGINAL winding for the retained face.
% - Quads are split along the better diagonal (shorter chord / better area).
% - N‑gons (>4) are fan‑triangulated as a best effort (warns once).
%
arguments
    meshes (1,:) cell
    names  (1,:) cell
    outFile (1,1) string
end
assert(numel(meshes)==numel(names), 'meshes and names must have same length.');

% ---------- small package parts ----------
contentTypes = ['<?xml version="1.0" encoding="UTF-8"?>' ...
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' ...
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' ...
    '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>' ...
    '</Types>'];
relsXML = ['<?xml version="1.0" encoding="UTF-8"?>' ...
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' ...
    '<Relationship Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel" Target="3D/3dmodel.model"/>' ...
    '</Relationships>'];

fos = java.io.FileOutputStream(char(outFile));
bos = java.io.BufferedOutputStream(fos, 1048576);  % 1 MiB buffer (was 1<<20)
zos = java.util.zip.ZipOutputStream(bos);
zos.setLevel(java.util.zip.Deflater.BEST_SPEED);    % or use NO_COMPRESSION
% % fallback:
% zos.setLevel(1);  % 1 = BEST_SPEED, 0 = NO_COMPRESSION, 9 = BEST_COMPRESSION

% Write small parts
writeEntryUTF8(zos, '[Content_Types].xml', contentTypes);
writeEntryUTF8(zos, '_rels/.rels',        relsXML);

% ---------- stream /3D/3dmodel.model ----------
entry = java.util.zip.ZipEntry('3D/3dmodel.model');
zos.putNextEntry(entry);
w = @(s) writeUTF8Chunk(zos, s);

w('<?xml version="1.0" encoding="UTF-8"?>');
w('<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">');
w('<resources>');

n = numel(meshes);
hasGeom = false(1,n);
warnedNgon = false;

for i = 1:n
    nm = xmlEscape(names{i});
    m  = meshes{i};

    if isempty(m) || ~isstruct(m) || ~isfield(m,'vertices') || ~isfield(m,'faces')
        warnPart(i, 'missing fields; emitting empty object');
        writeEmptyObject(w, i, nm);
        continue;
    end

    try
        V = double(m.vertices);
        F = double(m.faces);

        % ---- coerce faces to triangles WITHOUT destroying winding ----
        [F, warnedNgon] = coerceFacesToTriangles(V, F, warnedNgon);

        % ---- sanitize (drop NaNs/Inf/out-of-range/duplicate/degenerate/zero-area) ----
        [V,F] = sanitizeMesh(V,F);

        if isempty(V) || isempty(F)
            warnPart(i, 'no valid triangles after cleaning; emitting empty');
            writeEmptyObject(w, i, nm);
            continue;
        end

        hasGeom(i) = true;

        % ---- write object ----
        w(sprintf('<object id="%d" type="model" name="%s"><mesh><vertices>', i, nm));

        % vertices (batched)
        batch = 200000;
        nv = size(V,1);
        for sidx = 1:batch:nv
            eidx = min(sidx+batch-1, nv);
            buf = strings(eidx-sidx+1,1);
            k = 1;
            for vi = sidx:eidx
                buf(k) = sprintf('<vertex x="%.9f" y="%.9f" z="%.9f"/>', V(vi,1), V(vi,2), V(vi,3));
                k = k + 1;
            end
            w(strjoin(buf,''));
        end
        w('</vertices><triangles>');

        % triangles (0-based)
        nt = size(F,1);
        for sidx = 1:batch:nt
            eidx = min(sidx+batch-1, nt);
            buf = strings(eidx-sidx+1,1);
            k = 1;
            for fi = sidx:eidx
                v1 = F(fi,1)-1; v2 = F(fi,2)-1; v3 = F(fi,3)-1;
                buf(k) = sprintf('<triangle v1="%d" v2="%d" v3="%d"/>', v1, v2, v3);
                k = k + 1;
            end
            w(strjoin(buf,''));
        end

        w('</triangles></mesh></object>');

    catch ME
        warnPart(i, sprintf('exception during write (%s); emitting empty', ME.message));
        writeEmptyObject(w, i, nm);
    end
end

% ---- composite wrapper so Cura keeps parts together ----
ids = find(hasGeom);
wrapperId = n + 1;
w(sprintf('<object id="%d" type="model" name="assembly"><components>', wrapperId));
for id = ids
    w(sprintf('<component objectid="%d"/>', id));
end
w('</components></object>');

% build: only composite
w('</resources><build>');
if ~isempty(ids), w(sprintf('<item objectid="%d"/>', wrapperId)); end
w('</build></model>');

zos.closeEntry();
zos.close(); fos.close();
fprintf('3MF written: %s\n', char(outFile));
end

% ---------- helpers ----------

function writeEmptyObject(w, id, nm)
w(sprintf('<object id="%d" type="model" name="%s"><mesh><vertices/>', id, xmlEscape(nm)));
w('<triangles/></mesh></object>');
end

function warnPart(i, msg)
warning('export_3mf:part%d: %s', i, msg);
end

function [Ftri, warnedNgon] = coerceFacesToTriangles(V, F, warnedNgon)
% Acceptable inputs:
%  - Nx3 already → pass through
%  - 1x(3M) or (3M)x1 linearized → reshape to (M×3)
%  - Nx4 quads → split along best diagonal (shorter / better area)
%  - NxK with K>4 → best-effort fan around first index (warn once)

sz = size(F);
if numel(F)==0
    Ftri = zeros(0,3);
    return;
end

if sz(2)==3
    Ftri = F;
    return;
end

if sz(2)==1 || sz(1)==1
    L = numel(F);
    if mod(L,3)~=0
        error('Faces vector length %d not divisible by 3; cannot infer triangles.', L);
    end
    Ftri = reshape(F, 3, []).';
    return;
end

K = sz(2);
if K==4
    n = sz(1);
    Ftri = zeros(2*n,3, class(F));
    % choose diagonal per face using geometry
    for r = 1:n
        quad = F(r,:);
        p1 = V(quad(1),:); p2 = V(quad(2),:); p3 = V(quad(3),:); p4 = V(quad(4),:);
        d13 = sum((p1-p3).^2);
        d24 = sum((p2-p4).^2);
        % build two options and pick the one with larger minimal area
        optA = [quad([1 2 3]); quad([1 3 4])]; % diagonal 1-3
        optB = [quad([1 2 4]); quad([2 3 4])]; % diagonal 2-4
        Aareas = triArea(V,optA);
        Bareas = triArea(V,optB);
        scoreA = min(Aareas);
        scoreB = min(Bareas);
        if scoreB>scoreA || (scoreA==scoreB && d24<d13)
            Ftri(2*r-1:2*r,:) = optB;
        else
            Ftri(2*r-1:2*r,:) = optA;
        end
    end
    return;
end

% General K-gon fan triangulation (best effort)
if ~warnedNgon
    warning('export_3mf:ngonFan', 'Faces with K>4 detected; using fan triangulation as best effort.');
    warnedNgon = true;
end
n = sz(1);
Ftri = zeros((K-2)*n, 3, class(F));
row = 1;
for r = 1:n
    base = F(r,1);
    for j = 2:(K-1)
        Ftri(row,:) = [base, F(r,j), F(r,j+1)];
        row = row + 1;
    end
end
end

function [V,F] = sanitizeMesh(V,F)
% 1) basic shape checks
if size(V,2)~=3, error('vertices must be Nx3'); end
if size(F,2)~=3, error('faces must be Nx3 after coercion'); end

% 2) remove NaN/Inf vertices, remap faces
badV = any(~isfinite(V),2);
if any(badV)
    map = zeros(size(V,1),1);
    map(~badV) = 1:nnz(~badV);
    invalidF = badV(F(:,1)) | badV(F(:,2)) | badV(F(:,3));
    F = F(~invalidF,:);
    V = V(~badV,:);
    if isempty(F), V=[]; return; end
    F = map(F);
end

% 3) ensure integer & in-range indices
F = round(F);
nv = size(V,1);
maskRange = all(F>=1 & F<=nv,2);
F = F(maskRange,:);

% 4) drop exact duplicate triangles WITHOUT changing winding
if ~isempty(F)
    key = sort(F,2);                 % order‑independent key
    [~, ia] = unique(key, 'rows', 'stable');
    F = F(ia,:);                      % keep original winding of first
end

% 5) drop degenerate by index
deg = (F(:,1)==F(:,2)) | (F(:,2)==F(:,3)) | (F(:,1)==F(:,3));
F = F(~deg,:);
if isempty(F), V=[]; return; end

% 6) drop zero-area triangles (geometry)
A = triArea(V,F);
F = F(A>0,:);
if isempty(F), V=[]; end
end

function A = triArea(V,F)
v1 = V(F(:,2),:) - V(F(:,1),:);
v2 = V(F(:,3),:) - V(F(:,1),:);
cp = cross(v1,v2,2);
A = 0.5*sqrt(sum(cp.^2,2));
end

function writeEntryUTF8(zos, name, textChar)
entry = java.util.zip.ZipEntry(char(name));
zos.putNextEntry(entry);
writeUTF8Chunk(zos, textChar);
zos.closeEntry();
end

function writeUTF8Chunk(zos, textChar)
if isempty(textChar), return; end
chunkChars = 1e6;
L = numel(textChar);
for s = 1:chunkChars:L
    e = min(s+chunkChars-1, L);
    piece = textChar(s:e);
    jstr  = java.lang.String(piece);
    bytes = jstr.getBytes('UTF-8');
    zos.write(bytes, 0, numel(bytes));
end
end

function s = xmlEscape(in)
s = regexprep(char(string(in)), {'&','<','>','"',''''}, ...
                               {'&amp;','&lt;','&gt;','&quot;','&apos;'});
end
