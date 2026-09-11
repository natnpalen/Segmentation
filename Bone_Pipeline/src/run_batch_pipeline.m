function results = run_batch_pipeline(rootFolder, varargin)
% RUN_BATCH_PIPELINE  Segment every DICOM series found under a root folder.
%
%   results = run_batch_pipeline(rootFolder)
%   results = run_batch_pipeline(rootFolder, 'Option', value, ...)
%
% Searches rootFolder for subfolders containing DICOM images, then runs the
% bone pipeline on each one in segmentation-only mode: bone masks, NIfTI
% volumes and STL meshes, with no cortical/cancellous split, no specimen
% packing and no figures.
%
% Expected layout (nested image folders are found automatically):
%
%   root/
%     156L-1/DICOMOBJ/...    -> case '156L-1'
%     156R-2/DICOMOBJ/...    -> case '156R-2'
%
% Outputs are sorted into folders by file type, with HU volumes kept apart
% from mask volumes — see 'Organize' below. A run that dies on one scan keeps
% going and records the failure; re-running skips cases that already
% completed unless 'Overwrite' is true.
%
% Name-value options
%   'OutputRoot'   : ''   output root ('' = <rootFolder>/bone_pipeline_batch)
%   'Organize'     : 'type'  how outputs are foldered:
%                     'type' - pooled across cases by file type:
%                              <OutputRoot>/nifti_mask/<case>_bone_01_mask.nii.gz
%                              <OutputRoot>/nifti_hu/<case>_bone_01_hu.nii.gz
%                              <OutputRoot>/stl_smooth/, stl_voxelized/, summaries/
%                     'case' - one folder per scan, split by type inside:
%                              <OutputRoot>/<case>/nifti_mask/bone_01_mask.nii.gz
%                     'flat' - everything for a scan in <OutputRoot>/<case>/
%   'MaxBones'     : 1    bones to keep per scan (1 for single-bone scans;
%                          [] keeps everything separate_bones finds)
%   'MinFiles'     : 5    minimum DICOM files for a folder to count as a series
%   'Include'      : ''   regexp; only run cases whose name matches
%   'Exclude'      : ''   regexp; skip cases whose name matches
%   'Overwrite'    : false  re-run cases that already have outputs
%   'DryRun'       : false  list the discovered series and stop
%   'SaveMat'      : false  also write the large pipeline_results.mat per case
%   'PipelineArgs' : {}   extra name-value pairs passed to run_bone_pipeline
%                          (e.g. {'MinBoneVolMM3', 300, 'TargetIsoMM', 0.5})
%
% Output
%   results : struct array, one entry per case, with fields name, path,
%             status ('ok' | 'failed' | 'skipped' | 'no_bones'), n_bones,
%             n_bones_found (before the MaxBones cap), n_markers,
%             volumes_mm3, mean_hu, preset (which segmentation pass won),
%             flags (per-bone quality flags), elapsed_s, outputDir, message.
%
% batch_summary.csv carries the pass and quality flags per bone plus a
% 'review' column, so scans the segmentation was unsure about can be pulled
% out and checked by eye instead of being trusted silently.

o = struct( ...
    'OutputRoot',   '', ...
    'Organize',     'type', ...
    'MaxBones',     1, ...
    'MinFiles',     5, ...
    'Include',      '', ...
    'Exclude',      '', ...
    'Overwrite',    true, ...
    'DryRun',       false, ...
    'SaveMat',      false, ...
    'PipelineArgs', {{}} ...
);
o = utils.parse_opts(o, varargin{:});

if ~isfolder(rootFolder)
    error('Root folder not found: %s', rootFolder);
end

% ---- Make sure the pipeline is on the path ----
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

if isempty(o.OutputRoot)
    outputRoot = fullfile(rootFolder, 'bone_pipeline_batch');
else
    outputRoot = o.OutputRoot;
end

organize = lower(char(o.Organize));
if ~ismember(organize, {'type', 'case', 'flat'})
    error('Organize must be ''type'', ''case'' or ''flat'' (got ''%s'').', organize);
end
% Cases run into a staging folder first, then get filed by type, so a case
% folder never collides with a type folder.
stagingRoot = fullfile(outputRoot, '_staging');

fprintf('\n');
fprintf('==========================================================\n');
fprintf('  BONE SEGMENTATION — BATCH MODE\n');
fprintf('==========================================================\n');
fprintf('  Root   : %s\n', rootFolder);
fprintf('  Output : %s\n', outputRoot);
fprintf('  Layout : %s\n', layout_description(organize));
fprintf('==========================================================\n\n');

% ---- Discover DICOM series ----
fprintf('Scanning for DICOM series...\n');
t_scan = tic;
series = dicom.find_series_dirs(rootFolder, ...
    'MinFiles', o.MinFiles, 'Exclude', {outputRoot});
fprintf('Found %d DICOM series (%.1fs)\n\n', numel(series), toc(t_scan));

if isempty(series)
    warning(['No DICOM series found under %s. Check the path, or lower ' ...
             '''MinFiles'' if the scans have very few slices.'], rootFolder);
    results = empty_results();
    return;
end

% ---- Apply include/exclude filters ----
keep = true(1, numel(series));
for i = 1:numel(series)
    if ~isempty(o.Include) && isempty(regexpi(series(i).name, o.Include, 'once'))
        keep(i) = false;
    end
    if ~isempty(o.Exclude) && ~isempty(regexpi(series(i).name, o.Exclude, 'once'))
        keep(i) = false;
    end
end
if ~all(keep)
    fprintf('Filters kept %d of %d series\n\n', sum(keep), numel(series));
    series = series(keep);
end

for i = 1:numel(series)
    fprintf('  %3d. %-28s  %5d files  %s\n', i, series(i).name, ...
        series(i).n_files, series(i).relpath);
end
fprintf('\n');

if o.DryRun
    fprintf('Dry run — nothing processed.\n\n');
    results = empty_results();
    for i = 1:numel(series)
        results(i) = make_result(series(i), 'skipped', 'dry run'); %#ok<AGROW>
    end
    return;
end

if ~exist(outputRoot, 'dir'), mkdir(outputRoot); end

% ---- Run every case ----
n_cases = numel(series);
results = empty_results();
t_batch = tic;

for i = 1:n_cases
    s = series(i);

    if strcmp(organize, 'flat')
        caseOut = fullfile(outputRoot, s.name);   % pipeline writes final files
        finalOut = caseOut;
    else
        caseOut = fullfile(stagingRoot, s.name);  % pipeline writes, then filed
        if strcmp(organize, 'case')
            finalOut = fullfile(outputRoot, s.name);
        else
            finalOut = outputRoot;
        end
    end

    fprintf('----------------------------------------------------------\n');
    fprintf('[%d/%d] %s\n', i, n_cases, s.name);
    fprintf('----------------------------------------------------------\n');

    if ~o.Overwrite && case_is_done(outputRoot, s.name, organize)
        fprintf('  Already processed — skipping (use ''Overwrite'', true to redo)\n\n');
        r = make_result(s, 'skipped', 'existing outputs');
        r.outputDir = finalOut;
        results(i) = r; %#ok<AGROW>
        continue;
    end

    % Clear any partial staging left by an interrupted run
    if ~strcmp(organize, 'flat') && isfolder(caseOut)
        [ok, msg] = rmdir(caseOut, 's');
        if ~ok
            warning('Could not clear staging folder %s: %s', caseOut, msg);
        end
    end

    t_case = tic;
    try
        out = run_bone_pipeline(s.path, '', ...
            'CorticalCancellous', false, ...
            'PackSpecimens',      false, ...
            'ShowViewer',         false, ...
            'SaveOutputs',        true, ...
            'SaveMat',            o.SaveMat, ...
            'MaxBones',           o.MaxBones, ...
            'OutputDir',          caseOut, ...
            o.PipelineArgs{:});

        % File the outputs by type
        if ~strcmp(organize, 'flat')
            n_moved = utils.organize_outputs(caseOut, outputRoot, s.name, organize);
            fprintf('       Filed %d output files by type\n\n', n_moved);
        end

        r = make_result(s, 'ok', '');
        r.outputDir = finalOut;
        r.elapsed_s = toc(t_case);
        r.n_bones = numel(out.separation.bones);
        r.n_markers = out.separation.n_tags;
        if isfield(out.separation, 'n_bones_found')
            r.n_bones_found = out.separation.n_bones_found;
        else
            r.n_bones_found = r.n_bones;
        end
        if r.n_bones_found > r.n_bones
            r.message = sprintf('%d objects found, kept %d largest', ...
                r.n_bones_found, r.n_bones);
        end
        if isfield(out.separation, 'preset')
            r.preset = out.separation.preset;
        end
        if r.n_bones > 0
            r.volumes_mm3 = cellfun(@(b) b.volume_mm3, out.separation.bones);
            r.mean_hu = cellfun(@(b) b.mean_hu, out.separation.bones);
            r.flags = cell(1, r.n_bones);
            for bi = 1:r.n_bones
                if isfield(out.separation.bones{bi}, 'quality')
                    r.flags{bi} = strjoin(out.separation.bones{bi}.quality.flags, '|');
                else
                    r.flags{bi} = '';
                end
            end
        else
            r.status = 'no_bones';
            r.message = 'no bones found';
        end
        results(i) = r; %#ok<AGROW>

    catch ME
        fprintf(2, '  FAILED: %s\n\n', ME.message);
        r = make_result(s, 'failed', ME.message);
        r.outputDir = caseOut;   % partial files stay in staging for inspection
        r.elapsed_s = toc(t_case);
        results(i) = r; %#ok<AGROW>
    end

    % Write the summary after every case so a long run stays inspectable
    % (and survives an interrupted session).
    try
        write_batch_summary(outputRoot, rootFolder, results, o);
    catch ME
        warning('Batch summary write failed: %s', ME.message);
    end
end

elapsed = toc(t_batch);

% ---- Drop the staging folder if every case filed cleanly ----
if ~strcmp(organize, 'flat') && isfolder(stagingRoot)
    leftover = dir(stagingRoot);
    leftover = leftover(~ismember({leftover.name}, {'.', '..'}));
    if isempty(leftover)
        rmdir(stagingRoot);
    else
        fprintf('\n  Note: partial outputs from failed scans left in %s\n', stagingRoot);
    end
end

% ---- Final report ----
n_ok      = sum(strcmp({results.status}, 'ok'));
n_failed  = sum(strcmp({results.status}, 'failed'));
n_skipped = sum(strcmp({results.status}, 'skipped'));
n_none    = sum(strcmp({results.status}, 'no_bones'));

fprintf('==========================================================\n');
fprintf('  BATCH COMPLETE  (%.1f min)\n', elapsed / 60);
fprintf('==========================================================\n');
fprintf('  %d processed, %d skipped, %d with no bones, %d failed\n\n', ...
    n_ok, n_skipped, n_none, n_failed);

fprintf('  %-28s  %-9s  %6s  %11s  %8s\n', 'Case', 'Status', 'Bones', 'Volume', 'Time');
fprintf('  %-28s  %-9s  %6s  %11s  %8s\n', '----', '------', '-----', '------', '----');
for i = 1:numel(results)
    r = results(i);
    if isempty(r.volumes_mm3)
        vol_str = '-';
    else
        vol_str = sprintf('%.0f mm3', sum(r.volumes_mm3));
    end
    fprintf('  %-28s  %-9s  %6d  %11s  %6.1fs\n', ...
        trunc(r.name, 28), r.status, r.n_bones, vol_str, r.elapsed_s);
end

if n_failed > 0
    fprintf('\n  Failures:\n');
    for i = 1:numel(results)
        if strcmp(results(i).status, 'failed')
            fprintf('    %-28s  %s\n', trunc(results(i).name, 28), results(i).message);
        end
    end
end

fprintf('\n  Output: %s\n', outputRoot);
fprintf('==========================================================\n\n');
end


% =========================================================================
%  Local helper functions
% =========================================================================

function results = empty_results()
    results = struct('name', {}, 'path', {}, 'relpath', {}, 'status', {}, ...
        'message', {}, 'n_bones', {}, 'n_bones_found', {}, 'n_markers', {}, ...
        'volumes_mm3', {}, 'mean_hu', {}, 'preset', {}, 'flags', {}, ...
        'elapsed_s', {}, 'outputDir', {});
end


function r = make_result(s, status, message)
    r = struct();
    r.name        = s.name;
    r.path        = s.path;
    r.relpath     = s.relpath;
    r.status      = status;
    r.message     = message;
    r.n_bones       = 0;
    r.n_bones_found = 0;
    r.n_markers     = 0;
    r.volumes_mm3 = [];
    r.mean_hu     = [];
    r.preset      = '';
    r.flags       = {};
    r.elapsed_s   = 0;
    r.outputDir   = '';
end


function tf = case_is_done(outputRoot, caseName, organize)
% A case counts as finished once its summary and at least one bone mask are
% in place, wherever the chosen layout puts them.
    switch organize
        case 'flat'
            caseDir = fullfile(outputRoot, caseName);
            summaryFile = fullfile(caseDir, 'pipeline_summary.txt');
            maskGlob = fullfile(caseDir, 'bone_*_mask.nii.gz');
        case 'case'
            caseDir = fullfile(outputRoot, caseName);
            summaryFile = fullfile(caseDir, 'summaries', 'pipeline_summary.txt');
            maskGlob = fullfile(caseDir, 'nifti_mask', 'bone_*_mask.nii.gz');
        otherwise  % 'type'
            summaryFile = fullfile(outputRoot, 'summaries', ...
                [caseName '_pipeline_summary.txt']);
            maskGlob = fullfile(outputRoot, 'nifti_mask', ...
                [caseName '_bone_*_mask.nii.gz']);
    end

    tf = exist(summaryFile, 'file') == 2 && ~isempty(dir(maskGlob));
end


function s = layout_description(organize)
    switch organize
        case 'flat'
            s = 'one folder per scan (files together)';
        case 'case'
            s = 'one folder per scan, split by file type';
        otherwise
            s = 'pooled by file type across all scans';
    end
end


function write_batch_summary(outputRoot, rootFolder, results, o)
% Text summary plus a CSV with one row per bone, for downstream analysis.

    fid = fopen(fullfile(outputRoot, 'batch_summary.txt'), 'w');
    if fid < 0
        error('Cannot write batch summary in %s', outputRoot);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'BONE SEGMENTATION BATCH SUMMARY\n');
    fprintf(fid, '===============================\n');
    fprintf(fid, 'Date     : %s\n', datestr(now));
    fprintf(fid, 'Root     : %s\n', rootFolder);
    fprintf(fid, 'Output   : %s\n', outputRoot);
    fprintf(fid, 'Layout   : %s (%s)\n', lower(char(o.Organize)), ...
        layout_description(lower(char(o.Organize))));
    fprintf(fid, 'MaxBones : %s\n\n', mat2str(o.MaxBones));

    for i = 1:numel(results)
        r = results(i);
        fprintf(fid, '%-28s  %-9s  %d bone(s)  %.1fs\n', ...
            r.name, r.status, r.n_bones, r.elapsed_s);
        fprintf(fid, '    source: %s\n', r.path);
        for bi = 1:numel(r.volumes_mm3)
            if bi <= numel(r.flags) && ~isempty(r.flags{bi})
                fl = sprintf(' [%s]', strrep(r.flags{bi}, '|', ', '));
            else
                fl = '';
            end
            fprintf(fid, '    bone %d: %.1f mm3, mean HU %.0f%s\n', ...
                bi, r.volumes_mm3(bi), r.mean_hu(bi), fl);
        end
        if ~isempty(r.message)
            fprintf(fid, '    note: %s\n', r.message);
        end
    end

    cfid = fopen(fullfile(outputRoot, 'batch_summary.csv'), 'w');
    if cfid < 0
        error('Cannot write batch CSV in %s', outputRoot);
    end
    ccleanup = onCleanup(@() fclose(cfid)); %#ok<NASGU>
    fprintf(cfid, ['case,status,bone_index,volume_mm3,mean_hu,n_bones_found,' ...
                   'n_markers,pass,quality_flags,review,elapsed_s,' ...
                   'dicom_folder,output_folder\n']);
    for i = 1:numel(results)
        r = results(i);
        if isempty(r.volumes_mm3)
            fprintf(cfid, '%s,%s,,,,%d,%d,%s,,,%.1f,%s,%s\n', csv(r.name), csv(r.status), ...
                r.n_bones_found, r.n_markers, csv(r.preset), r.elapsed_s, ...
                csv(r.path), csv(r.outputDir));
        else
            for bi = 1:numel(r.volumes_mm3)
                if bi <= numel(r.flags), fl = r.flags{bi}; else, fl = ''; end
                fprintf(cfid, '%s,%s,%d,%.1f,%.0f,%d,%d,%s,%s,%d,%.1f,%s,%s\n', ...
                    csv(r.name), csv(r.status), bi, r.volumes_mm3(bi), ...
                    r.mean_hu(bi), r.n_bones_found, r.n_markers, ...
                    csv(r.preset), csv(fl), ~isempty(fl), r.elapsed_s, ...
                    csv(r.path), csv(r.outputDir));
            end
        end
    end
end


function s = csv(s)
    s = strrep(char(s), '"', '""');
    if any(ismember(s, ',"'))
        s = sprintf('"%s"', s);
    end
end


function s = trunc(s, n)
    s = char(s);
    if numel(s) > n, s = ['...' s(end-n+4:end)]; end
end
