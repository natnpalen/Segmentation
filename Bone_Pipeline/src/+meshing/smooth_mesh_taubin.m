function V = smooth_mesh_taubin(V, F, iters, lambda, mu)
% SMOOTH_MESH_TAUBIN  Volume-preserving surface smoothing.
%
%   V = meshing.smooth_mesh_taubin(V, F)
%   V = meshing.smooth_mesh_taubin(V, F, iters, lambda, mu)
%
% Plain Laplacian smoothing pulls every vertex toward its neighbours, so the
% mesh shrinks a little on every pass — run it long enough and a bone loses
% real surface, ridges round off and the STL no longer matches the mask it
% came from. Taubin smoothing alternates the same shrinking step (lambda)
% with a slightly larger inflating step (mu, negative), which cancels the
% shrinkage while still removing voxel staircase.
%
% Inputs
%   V      : Nx3 vertices (mm)
%   F      : Mx3 faces, 1-based
%   iters  : number of lambda/mu pairs (default 8). 0 returns V unchanged.
%   lambda : shrink weight, 0..1 (default 0.5)
%   mu     : inflate weight, negative and slightly larger in magnitude than
%            lambda (default -0.53)
%
% Reference: Taubin, "A signal processing approach to fair surface design",
% SIGGRAPH 1995.

if nargin < 3 || isempty(iters),  iters = 8;      end
if nargin < 4 || isempty(lambda), lambda = 0.5;   end
if nargin < 5 || isempty(mu),     mu = -0.53;     end

if iters <= 0 || isempty(V) || isempty(F)
    return;
end

V = double(V);
F = double(F);
nv = size(V, 1);

% ---- Vertex adjacency (built in one shot, not face by face) ----
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = [E; E(:, [2 1])];                       % symmetric
adj = sparse(E(:,1), E(:,2), 1, nv, nv);
adj = spones(adj);                          % collapse duplicate edges

valence = full(sum(adj, 2));
valence(valence == 0) = 1;                  % isolated vertices stay put

for k = 1:iters
    % Shrinking pass
    delta = (adj * V) ./ valence - V;
    V = V + lambda * delta;

    % Inflating pass — undoes the volume loss of the pass above
    delta = (adj * V) ./ valence - V;
    V = V + mu * delta;
end
end
