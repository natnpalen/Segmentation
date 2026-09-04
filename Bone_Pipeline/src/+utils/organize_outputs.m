function moved = organize_outputs(srcDir, destRoot, caseName, mode)
% ORGANIZE_OUTPUTS  Sort one case's pipeline outputs into per-file-type folders.
%
%   moved = utils.organize_outputs(srcDir, destRoot, caseName, mode)
%
% Moves every file the pipeline wrote for one scan out of srcDir and into a
% folder named for its type, keeping HU volumes apart from mask volumes.
%
%   mode 'type' : destRoot/<type>/<caseName>_<file>
%                 Every case pools into the same type folders, so all masks
%                 (or all STLs) across the whole batch sit together. The
%                 case name is prefixed onto each file so it stays traceable.
%
%   mode 'case' : destRoot/<caseName>/<type>/<file>
%                 One folder per scan, split by type inside it.
%
% File types:
%   stl_smooth, stl_voxelized, stl, nifti_hu, nifti_mask, nifti_cortical,
%   nifti_cancellous, nifti_other, summaries, mat, figures, other
%
% srcDir is removed once it is empty. Returns the number of files moved; a
% missing or empty srcDir is not an error and returns 0.

moved = 0;
if ~isfolder(srcDir), return; end

mode = lower(char(mode));
if ~ismember(mode, {'type', 'case'})
    error('organize_outputs: mode must be ''type'' or ''case''.');
end

listing = dir(srcDir);
listing = listing(~ismember({listing.name}, {'.', '..'}));
files = listing(~[listing.isdir]);

for i = 1:numel(files)
    name = files(i).name;
    folder = type_folder(name);

    if strcmp(mode, 'type')
        destDir = fullfile(destRoot, folder);
        destName = [caseName '_' name];
    else
        destDir = fullfile(destRoot, caseName, folder);
        destName = name;
    end

    if ~exist(destDir, 'dir'), mkdir(destDir); end

    src = fullfile(srcDir, name);
    dst = fullfile(destDir, destName);
    [ok, msg] = movefile(src, dst, 'f');
    if ok
        moved = moved + 1;
    else
        warning('Could not move %s to %s: %s', src, dst, msg);
    end
end

% Drop the now-empty source folder (left alone if anything remains)
remaining = dir(srcDir);
remaining = remaining(~ismember({remaining.name}, {'.', '..'}));
if isempty(remaining)
    [ok, msg] = rmdir(srcDir);
    if ~ok
        warning('Could not remove %s: %s', srcDir, msg);
    end
end
end


% =========================================================================
%  FILE TYPE CLASSIFICATION
% =========================================================================
function f = type_folder(name)
    n = lower(name);

    if endsWith(n, '_smooth.stl')
        f = 'stl_smooth';
    elseif endsWith(n, '_voxelized.stl')
        f = 'stl_voxelized';
    elseif endsWith(n, '.stl')
        f = 'stl';
    elseif endsWith(n, {'_hu.nii.gz', '_hu.nii'})
        f = 'nifti_hu';
    elseif endsWith(n, {'_mask.nii.gz', '_mask.nii'})
        f = 'nifti_mask';
    elseif endsWith(n, {'_cortical.nii.gz', '_cortical.nii'})
        f = 'nifti_cortical';
    elseif endsWith(n, {'_cancellous.nii.gz', '_cancellous.nii'})
        f = 'nifti_cancellous';
    elseif endsWith(n, {'.nii.gz', '.nii'})
        f = 'nifti_other';
    elseif endsWith(n, '.mat')
        f = 'mat';
    elseif endsWith(n, {'.png', '.fig', '.jpg', '.tif'})
        f = 'figures';
    elseif endsWith(n, {'.txt', '.csv', '.json'})
        f = 'summaries';
    else
        f = 'other';
    end
end
