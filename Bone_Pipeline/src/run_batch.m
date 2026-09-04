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

clear all; close all; clc; %#ok<CLALL>

% ---- Set your paths here ----
rootFolder = 'C:\Users\natha\OneDrive\Documents\Nathaniel\For Nick\3D Scaphoid\Dicom Series Scaphoid';
outputRoot = '';   % '' = <rootFolder>\bone_pipeline_batch

% ---- Add this pipeline to the path ----
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ---- Run ----
% MaxBones = 1 keeps only the largest bone per scan (these scans hold one
% bone and one marker). Set it to [] to keep every bone that is found.
results = run_batch_pipeline(rootFolder, ...
    'OutputRoot', outputRoot, ...
    'MaxBones',   1);

% Tip: add 'DryRun', true to list what would be processed without running,
% or 'Overwrite', true to redo cases that already have outputs.
