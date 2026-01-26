function phi = gyroid_phi_trilinear(p, rhoGrid, origin, spacing, w, rho_eps)
% Evaluate gyroid φ at point p using trilinear interpolation of the PHI-FIELD,
% as described in Wegner & Campbell, 2024. This is the corrected approach.
%
%  1) Find the 8 corners surrounding point p.
%  2) For each corner, calculate its local cell size L from its density rho.
%  3) For each corner, evaluate phi(p, L_corner). This gives 8 distinct phi values.
%  4) Trilinearly interpolate these 8 phi values to get the final result at p.

sz = size(rhoGrid);

% Quick reject if p is outside the rhoGrid bbox (1/2 voxel tolerance)
minPt = origin;
maxPt = origin + spacing.*(sz - 1);
if any(p < minPt - 1e-9) || any(p > maxPt + 1e-9)
    phi = NaN; return;
end

% --- Step 1: Find 8 corners and interpolation weights ---
% This part remains the same as your original function.
u  = (p - origin)./spacing + 1;     % 1-based fractional index
i0 = floor(u);
i0 = max([1 1 1], min(sz - 1, i0)); % Clamp to ensure i0 and i1 are in-bounds
i1 = i0 + 1;
t  = u - i0;                        % Interpolation factors
tx = t(1); ty = t(2); tz = t(3);

% Fetch 8 corner densities
rho000 = rhoGrid(i0(1), i0(2), i0(3));
rho100 = rhoGrid(i1(1), i0(2), i0(3));
rho010 = rhoGrid(i0(1), i1(2), i0(3));
rho110 = rhoGrid(i1(1), i1(2), i0(3));
rho001 = rhoGrid(i0(1), i0(2), i1(3));
rho101 = rhoGrid(i1(1), i0(2), i1(3));
rho011 = rhoGrid(i0(1), i1(2), i1(3));
rho111 = rhoGrid(i1(1), i1(2), i1(3));
rho8 = [rho000 rho100 rho010 rho110 rho001 rho101 rho011 rho111];

% --- Step 2: Calculate L for all 8 corners and then evaluate φ at point p ---
% This is the core change from the original function.

% Guardrail the densities and convert all 8 to cell sizes (L)
rho8_guarded = max(rho8, rho_eps);
L8 = 2.432 * w ./ rho8_guarded;

% Evaluate phi at the target point 'p' using each of the 8 different L-values
phi8 = zeros(1, 8);
for i = 1:8
    % The point 'p' stays constant; only the cell size 'L' changes.
    phi8(i) = infill.gyroid_phi_at_point(p, L8(i));
end

% --- Step 3: Trilinearly interpolate the resulting 8 φ values ---
% We reuse the robust NaN-handling logic, but apply it to the final phi blend.
% The paper describes this interpolation of the function's value[cite: 91, 110].

% Standard trilinear weights (same as before)
w000 = (1-tx)*(1-ty)*(1-tz);
w100 =    tx *(1-ty)*(1-tz);
w010 = (1-tx)* ty *(1-tz);
w110 =    tx * ty *(1-tz);
w001 = (1-tx)*(1-ty)* tz;
w101 =    tx *(1-ty)* tz;
w011 = (1-tx)* ty * tz;
w111 =    tx * ty * tz;
w8   = [w000   w100   w010   w110   w001   w101   w011   w111];

% NaN-robust blend: use the validity of the original rho values to guide the blend
valid = isfinite(rho8);
if ~any(valid)
    phi = NaN; return;                 % All surrounding corners are outside the mask
end

w8(~valid)   = 0;   % Ignore weights for corners with NaN density
phi8(~valid) = 0;   % Set corresponding phi to 0 so it doesn't contribute NaN

wsum = sum(w8);
if wsum <= 0
    % This can happen if the only valid corners have zero weight
    phi = NaN; return;
end

% The final value is the weighted average of the 8 phi results
phi = sum(w8 .* phi8) / wsum;

end