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
% Each case is written to <OutputRoot>/<caseName>/. A run that dies on one
% scan keeps going and records the failure; re-running skips cases that
% already completed unless 'Overwrite' is true.
%
% Name-value options
%   'OutputRoot'   : ''   output root ('' = <rootFolder>/bone_pipeline_batch)
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
%             volumes_mm3, mean_hu, elapsed_s, outputDir, message.

o = struct( ...
    'OutputRoot',   '', ...
    'MaxBones',     1, ...
    'MinFiles',     5, ...
    'Include',      '', ...
    'Exclude',      '', ...
    'Overwrite',    false, ...
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

fprintf('\n');
fprintf('==========================================================\n');
fprintf('  BONE SEGMENTATION — BATCH MODE\n');
fprintf('==========================================================\n');
fprintf('  Root   : %s\n', rootFolder);
fprintf('  Output : %s\n', outputRoot);
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
    caseOut = fullfile(outputRoot, s.name);

    fprintf('----------------------------------------------------------\n');
    fprintf('[%d/%d] %s\n', i, n_cases, s.name);
    fprintf('----------------------------------------------------------\n');

    if ~o.Overwrite && case_is_done(caseOut)
        fprintf('  Already processed — skipping (use ''Overwrite'', true to redo)\n\n');
        r = make_result(s, 'skipped', 'existing outputs');
        r.outputDir = caseOut;
        results(i) = r; %#ok<AGROW>
        continue;
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

        r = make_result(s, 'ok', '');
        r.outputDir = caseOut;
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
        if r.n_bones > 0
            r.volumes_mm3 = cellfun(@(b) b.volume_mm3, out.separation.bones);
            r.mean_hu = cellfun(@(b) b.mean_hu, out.separation.bones);
        else
            r.status = 'no_bones';
            r.message = 'no bones found';
        end
        results(i) = r; %#ok<AGROW>

    catch ME
        fprintf(2, '  FAILED: %s\n\n', ME.message);
        r = make_result(s, 'failed', ME.message);
        r.outputDir = caseOut;
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
        'volumes_mm3', {}, 'mean_hu', {}, 'elapsed_s', {}, 'outputDir', {});
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
    r.elapsed_s   = 0;
    r.outputDir   = '';
end


function tf = case_is_done(caseOut)
% A case counts as finished once the pipeline wrote its summary and at
% least one bone mask.
    tf = isfolder(caseOut) && ...
         exist(fullfile(caseOut, 'pipeline_summary.txt'), 'file') == 2 && ...
         ~isempty(dir(fullfile(caseOut, 'bone_*_mask.nii.gz')));
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
    fprintf(fid, 'MaxBones : %s\n\n', mat2str(o.MaxBones));

    for i = 1:numel(results)
        r = results(i);
        fprintf(fid, '%-28s  %-9s  %d bone(s)  %.1fs\n', ...
            r.name, r.status, r.n_bones, r.elapsed_s);
        fprintf(fid, '    source: %s\n', r.path);
        for bi = 1:numel(r.volumes_mm3)
            fprintf(fid, '    bone %d: %.1f mm3, mean HU %.0f\n', ...
                bi, r.volumes_mm3(bi), r.mean_hu(bi));
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
                   'n_markers,elapsed_s,dicom_folder,output_folder\n']);
    for i = 1:numel(results)
        r = results(i);
        if isempty(r.volumes_mm3)
            fprintf(cfid, '%s,%s,,,,%d,%d,%.1f,%s,%s\n', csv(r.name), csv(r.status), ...
                r.n_bones_found, r.n_markers, r.elapsed_s, csv(r.path), csv(r.outputDir));
        else
            for bi = 1:numel(r.volumes_mm3)
                fprintf(cfid, '%s,%s,%d,%.1f,%.0f,%d,%d,%.1f,%s,%s\n', ...
                    csv(r.name), csv(r.status), bi, r.volumes_mm3(bi), ...
                    r.mean_hu(bi), r.n_bones_found, r.n_markers, r.elapsed_s, ...
                    csv(r.path), csv(r.outputDir));
            end
        end
    end
end


function s = csv(s)
    s = strrep(char(s), '"', '""');
    if any(ismember(s, ',"'))
        s = ['"' s '"'];
    end
end


function s = trunc(s, n)
    s = char(s);
    if numel(s) > n, s = ['...' s(end-n+4:end)]; end
end
