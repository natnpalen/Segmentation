% RUN_BATCH  Quick-start script for batch bone segmentation.
%
% Point rootFolder at the folder that holds one subfolder per scan, then hit
% Run (F5). Every DICOM series underneath it is segmented and exported as
% NIfTI masks + STL meshes — no cortical/cancellous split, no specimen
% packing, no figures.
%
%   root/
%     156L-1/DICOMOBJ/...
%     156R-2/DICOMOBJ/...
%     ...
%
% Results are filed by file type, HU volumes kept apart from masks:
%
%   bone_pipeline_batch/
%     nifti_mask/      156L-1_bone_01_mask.nii.gz, 156R-2_...
%     nifti_hu/        156L-1_bone_01_hu.nii.gz, ...
%     stl_smooth/      156L-1_bone_01_smooth.stl, ...
%     stl_voxelized/   156L-1_bone_01_voxelized.stl, ...
%     summaries/       156L-1_pipeline_summary.txt, ...
%     batch_summary.txt / batch_summary.csv

clear all; close all; clc; %#ok<CLALL>

% ---- Set your paths here ----
rootFolder = 'C:\Users\natha\OneDrive\Documents\Nathaniel\For Nick\New Bone Scans';
outputRoot = '';   % '' = <rootFolder>\bone_pipeline_batch

% ---- Add this pipeline to the path ----
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ---- Run ----
% MaxBones = 1 keeps only the largest bone per scan (these scans hold one
% bone and one marker). Set it to [] to keep every bone that is found.
results = run_batch_pipeline(rootFolder, ...
    'OutputRoot', outputRoot, ...
    'MaxBones',   1, ...
    'Organize',   'type');

% Tip: add 'DryRun', true to list what would be processed without running,
% or 'Overwrite', true to redo cases that already have outputs.
%
% 'Organize' controls the output layout:
%   'type' - pooled by file type across all scans (above)
%   'case' - one folder per scan, split by file type inside it
%   'flat' - one folder per scan, all its files together
