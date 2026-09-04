function opts = segment_options(opts, preset)
% SEGMENT_OPTIONS  Density and tissue knobs for bone separation.
%
%   opts = bone.segment_options(opts)                % fill in defaults
%   opts = bone.segment_options(opts, 'lowdensity')  % apply a preset
%
% Bone separation is threshold-driven, and the thresholds are tuned for
% healthy cortical bone. These five knobs scale that tuning so a scan the
% defaults cannot handle can be retried under different assumptions.
%
% Knobs
%   DensityScale  : 1.0  multiplies every HU floor in core selection, the
%                        FMM sweep and boundary refinement. Below 1 keeps
%                        low-density (osteoporotic) bone that the default
%                        floors would carve away; above 1 tightens.
%   CorePrctile   : 94   percentile of local HU used for the dense-core
%                        seed. Lower it when the bone has no dense core.
%   FMMThreshMax  : 0.42 upper end of the front-arrival sweep. Higher lets
%                        the region grow further before scoring.
%   TissueScrub   : 1.0  multiplies the surface tissue-removal threshold.
%                        Above 1 strips more clinging soft tissue.
%   MinBoneHU     : 50   a candidate whose mean HU is below this is not
%                        bone and is discarded.
%
% Presets
%   'standard'    : the tuned defaults, unchanged. Every knob is neutral,
%                   so a standard pass is bit-identical to the pipeline
%                   before presets existed.
%   'lowdensity'  : relaxed floors for osteoporotic or low-contrast scans,
%                   where the standard pass finds nothing or only a
%                   fragment of the bone.
%   'tissue'      : tighter floors and a harder surface scrub, for scans
%                   where the standard pass swallowed soft tissue.
%
% Presets only change these knobs plus MinBoneVolMM3. Marker handling,
% seeding and geometry cleanup are the same in every pass.

DEFAULTS = struct( ...
    'DensityScale', 1.0, ...
    'CorePrctile',  94, ...
    'FMMThreshMax', 0.42, ...
    'TissueScrub',  1.0, ...
    'MinBoneHU',    50 ...
);

names = fieldnames(DEFAULTS);
for i = 1:numel(names)
    if ~isfield(opts, names{i}) || isempty(opts.(names{i}))
        opts.(names{i}) = DEFAULTS.(names{i});
    end
end

if nargin < 2 || isempty(preset), return; end

switch lower(char(preset))
    case 'standard'
        % Defaults as supplied — nothing to change.

    case 'lowdensity'
        opts.DensityScale = 0.45;
        opts.CorePrctile  = 85;
        opts.FMMThreshMax = 0.55;
        opts.TissueScrub  = 1.0;
        opts.MinBoneHU    = 35;
        if isfield(opts, 'MinBoneVolMM3')
            opts.MinBoneVolMM3 = opts.MinBoneVolMM3 * 0.5;
        end

    case 'tissue'
        opts.DensityScale = 1.25;
        opts.CorePrctile  = 96;
        opts.FMMThreshMax = 0.35;
        opts.TissueScrub  = 1.7;
        opts.MinBoneHU    = 80;

    otherwise
        error('Unknown segmentation preset: %s', preset);
end
end
