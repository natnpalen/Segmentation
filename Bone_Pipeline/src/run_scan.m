% RUN_SCAN  Quick-start script for the bone segmentation pipeline.
%
% Set the DICOM and STL folder paths below, then hit Run (F5).

clear all; close all; clc; %#ok<CLALL>

% ---- Set your paths here ----
dicomFolder = 'C:\Users\natha\Downloads\Metacarpals from 8-5-26-20260814T194951Z-1-001\Metacarpals from 8-5-26\Batch 1\ScalarVolume_58';
stlFolder   = 'C:\Users\natha\OneDrive\Documents\Nathaniel\For Nick\Mechancial Specimens';

% ---- Add this pipeline to the path ----
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ---- Run ----
% Optional modes (add as name-value pairs below):
%   'SplitCorticalCancellous', false  -> turn OFF the cortical/cancellous
%       sectioning; each bone is kept as one whole region and specimen
%       packing automatically runs whole-bone.
%   'ShavedBoneMode', true  -> for machined specimens, e.g. metacarpals
%       with cortical bone shaved flat for 3-point bending. Lowers the
%       minimum bone size (150 mm^3 instead of 500) and uses gentler
%       surface cleanup so thin cortical plates survive.
%   'PackingOrientations', 24  -> finer rotation sweep for specimen
%       packing (default 6 = 90-degree steps only). Helps elongated
%       specimens (Bend, Shear) fit inside curved bones.
out = run_bone_pipeline(dicomFolder, stlFolder, ...
    'PackSpecimens', false, ...
    'PackWholeBone', false, ....
    'ShavedBoneMode', true, ...
    'SplitCorticalCancellous', false, ...
    'SaveOutputs',   true);