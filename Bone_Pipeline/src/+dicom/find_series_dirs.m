function series = find_series_dirs(rootFolder, varargin)
% FIND_SERIES_DIRS  Locate DICOM series folders beneath a root folder.
%
%   series = dicom.find_series_dirs(rootFolder)
%   series = dicom.find_series_dirs(rootFolder, 'MinFiles', 10)
%
% Walks the folder tree under rootFolder and returns every directory that
% directly contains at least MinFiles DICOM images. Handles the common
% scanner-export layout where each specimen is its own subfolder holding a
% nested image folder:
%
%   root/
%     156L-1/DICOMOBJ/0000004F, 00000050, ...
%     156R-2/DICOMOBJ/...
%
% Files are identified by the 'DICM' magic bytes at offset 128, so
% extensionless scanner exports are recognized. Folders whose files lack the
% preamble fall back to a dicominfo probe on a few files.
%
% Name-value options
%   'MinFiles' : 5    minimum DICOM files for a folder to count as a series
%   'MaxDepth' : 8    maximum folder depth to descend below rootFolder
%   'Exclude'  : {}   folder paths to skip (e.g. the batch output folder)
%
% Output
%   series : struct array with fields
%     .path     full path to the series folder
%     .relpath  path relative to rootFolder
%     .name     suggested case name (unique across the returned array)
%     .n_files  number of DICOM files found directly in the folder

o = struct('MinFiles', 5, 'MaxDepth', 8, 'Exclude', {{}});
o = utils.parse_opts(o, varargin{:});

if ~isfolder(rootFolder)
    error('Root folder not found: %s', rootFolder);
end
rootFolder = char(canonical_path(rootFolder));

excl = o.Exclude;
if ischar(excl) || isstring(excl), excl = {char(excl)}; end
for i = 1:numel(excl)
    excl{i} = char(canonical_path(excl{i}));
end

series = struct('path', {}, 'relpath', {}, 'name', {}, 'n_files', {});
series = scan_dir(rootFolder, rootFolder, 0, o, excl, series);

if isempty(series), return; end

% Sort by relative path so batch runs are reproducible
[~, ord] = sort(lower({series.relpath}));
series = series(ord);

% Assign unique, human-meaningful case names
names = cell(1, numel(series));
for i = 1:numel(series)
    names{i} = case_name(series(i).relpath);
end
names = make_unique(names);
for i = 1:numel(series)
    series(i).name = names{i};
end
end


% =========================================================================
%  RECURSIVE WALK
% =========================================================================
function series = scan_dir(folder, rootFolder, depth, o, excl, series)
    for i = 1:numel(excl)
        if strcmpi(folder, excl{i}) || starts_with_path(folder, excl{i})
            return;
        end
    end

    % Never descend into this pipeline's own results (unless the caller
    % pointed the search straight at one)
    if depth > 0
        [~, leaf] = fileparts(folder);
        if any(strcmpi(leaf, {'bone_pipeline_batch', 'bone_pipeline_outputs'}))
            return;
        end
    end

    listing = dir(folder);
    listing = listing(~ismember({listing.name}, {'.', '..'}));

    files = listing(~[listing.isdir]);
    n_dicom = count_dicom_files(folder, files, o.MinFiles);

    if n_dicom >= o.MinFiles
        rec = struct();
        rec.path = folder;
        rec.relpath = relative_path(folder, rootFolder);
        rec.name = '';
        rec.n_files = n_dicom;
        series(end+1) = rec; %#ok<AGROW>
    end

    if depth >= o.MaxDepth, return; end

    subdirs = listing([listing.isdir]);
    for k = 1:numel(subdirs)
        series = scan_dir(fullfile(folder, subdirs(k).name), rootFolder, ...
            depth + 1, o, excl, series);
    end
end


% =========================================================================
%  DICOM FILE COUNTING
% =========================================================================
function n = count_dicom_files(folder, files, min_files)
    n = 0;
    if isempty(files), return; end

    for i = 1:numel(files)
        if files(i).bytes < 512, continue; end
        if has_dicm_magic(fullfile(folder, files(i).name))
            n = n + 1;
        end
    end

    if n > 0, return; end

    % Fallback: DICOM files written without the 128-byte preamble. Probe a
    % handful with dicominfo; if they read, count every plausible file.
    candidates = files([files.bytes] >= 512);
    if numel(candidates) < min_files, return; end

    n_probe = min(3, numel(candidates));
    ok = 0;
    for i = 1:n_probe
        try
            dicominfo(fullfile(folder, candidates(i).name));
            ok = ok + 1;
        catch
        end
    end
    if ok == n_probe
        n = numel(candidates);
    end
end


function tf = has_dicm_magic(filepath)
    tf = false;
    fid = fopen(filepath, 'r', 'ieee-le');
    if fid < 0, return; end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    if fseek(fid, 128, 'bof') ~= 0, return; end
    magic = fread(fid, 4, '*char');
    tf = numel(magic) == 4 && strcmp(magic(:).', 'DICM');
end


% =========================================================================
%  NAMING
% =========================================================================
function name = case_name(relpath)
% Build a readable case name from the relative path, dropping generic
% wrapper folder names so root/156L-1/DICOMOBJ becomes '156L-1'.

    GENERIC = {'dicom', 'dicomobj', 'dicomdat', 'dicomdir', 'images', ...
               'image', 'img', 'data', 'files', 'export'};

    % regexp split rather than strsplit: strsplit runs escape processing on
    % its delimiters, which mangles the Windows separator.
    parts = regexp(relpath, '[\\/]+', 'split');
    parts = parts(~cellfun(@isempty, parts));
    if isempty(parts), parts = {'series'}; end

    keep = parts(~ismember(lower(parts), GENERIC));
    if isempty(keep), keep = parts(end); end

    name = strjoin(keep, '_');
    name = regexprep(name, '[^\w\-\.]', '_');
    name = regexprep(name, '_+', '_');
    name = regexprep(name, '^_|_$', '');
    if isempty(name), name = 'series'; end
end


function names = make_unique(names)
    seen = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:numel(names)
        key = lower(names{i});
        if isKey(seen, key)
            seen(key) = seen(key) + 1;
            names{i} = sprintf('%s_%d', names{i}, seen(key));
        else
            seen(key) = 1;
        end
    end
end


% =========================================================================
%  PATH UTILITIES
% =========================================================================
function p = canonical_path(p)
    p = char(p);
    d = dir(p);
    if ~isempty(d) && isfield(d, 'folder') && isfolder(p)
        % dir('.') reports the resolved absolute folder of the first entry
        p = d(1).folder;
    end
    p = regexprep(p, '[\\/]+$', '');
end


function rel = relative_path(p, rootFolder)
    p = char(p);
    rootFolder = char(rootFolder);
    if strncmpi(p, rootFolder, numel(rootFolder))
        rel = p(numel(rootFolder)+1:end);
    else
        rel = p;
    end
    rel = regexprep(rel, '^[\\/]+', '');
    if isempty(rel), [~, leaf] = fileparts(rootFolder); rel = leaf; end
end


function tf = starts_with_path(p, prefix)
    prefix = [regexprep(prefix, '[\\/]+$', '') filesep];
    tf = strncmpi(p, prefix, numel(prefix));
end
