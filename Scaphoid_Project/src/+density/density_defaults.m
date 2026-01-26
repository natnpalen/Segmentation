function D = density_defaults()
% DENSITY_DEFAULTS  Centralized knobs for density analysis (used everywhere).

D = struct();
D.k               = 6;           % K-means zones
D.island_min_mm3  = 27.0;         % remove tiny islands (< this volume)
D.majority_iters  = 3;           % passes of neighborhood vote to fill gaps
D.neigh_mm        = 2.0;         % label smoothing radius (mm)
D.sigma_mm        = 0.5;         % masked Gaussian blur on HU (mm)
D.clamp_prct      = [2 98];      % robust HU clamp percentiles
D.nbrhood_size_vox= [6 6 6];     % majority voting window (voxels)
D.rngSeed         = [];          % [] → default RNG; or set numeric for reproducibility
end
