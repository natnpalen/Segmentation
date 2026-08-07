% RUN_SCAN  Quick-start script for the bone segmentation pipeline.
%
% Set the DICOM and STL folder paths below, then hit Run (F5).

clear all; close all; clc; %#ok<CLALL>

% ---- Set your paths here ----
dicomFolder = 'C:\Users\natha\OneDrive\Documents\Nathaniel\For Nick\New Bone Scans\156L-1\DICOMOBJ';
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
out = run_bone_pipeline(dicomFolder, stlFolder, ...
    'PackSpecimens', true, ...
    'PackWholeBone', true, ...
    'SaveOutputs',   true);
