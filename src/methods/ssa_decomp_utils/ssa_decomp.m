function modes = ssa_decomp(x, L, nc, epsilon, force_supk)
% function modes = ssa_decomp(x, L, nc, epsilon, force_supk)
%
% SSA with hierarchical clustering  for components grouping
% In the hierarchical clustering, the last components whose energy contribution is relatively
% null are eliminated from the clustering procedure. epsilon serves as threshold. 
%
% Input:
% x: input signal
% L: SSA window length
% nc: number of components to extract
% epsilon: threshold applied on the values in the singular spectrum (default: 3e-2)
% force_supk: constraint the number of singular values to use (default: supk depends on epsilon)
%
%
% Authors: J. Harmouche and D. Fourer
% Date: Apr. 2016
%
% Ref: [J. Harmouche, D.Fourer, F. Auger, P.Flandrin and P. Borgnat. The Sliding Singular Spectrum Analysis:
% a Data-Driven Non-Stationary Signal Decomposition Tool. IEEE TSP 2016.]


if ~exist('epsilon', 'var')
 epsilon = 3e-3; %0.5; %
end

if ~exist('force_supk', 'var')
 force_supk = inf;    
end


Nx=length(x);

%modes = zeros(Nx, 1);
if nc < 1
 %modes = zeros(Nx, 1);
 return;
end

[~,d]=ssa(x,L,1);

%supk = find( (d/max(d))<epsilon, 1) % C'est pour �liminer les composantes negligeables en amplitude, pouvant alterer le clustering 
supk = length(find( (d/max(d))>epsilon));


if supk <  nc
 warning('epsilon is too high, trying with eps=%.2f', d(nc)/max(d))
 supk = nc;
end

if force_supk < inf
 %fprintf(1, 'Forcing supk to %d', force_supk);
 supk = force_supk;
end

y = zeros(Nx,supk);
for I = 1:supk;
 [y(:,I),~]=ssa(x,L,I); 
end

%wcorr = pdist(y.', 'Euclidean');
%wcorr = 1-squareform(pdist(y.', 'correlation'));
wcorr = wCorrMat2(x,L,supk);

if ~isempty(wcorr)
  Z = linkage(wcorr);                       
  T = cluster(Z,'maxclust', nc); 
 %[~,~,T] = kmeans(y.', nc);  %, 'distance', 'correlation'
 
 modes = zeros(Nx, nc);

 for i = 1:nc       
  idx = find(T == i);

  if length(idx) > 1
   modes(:, i) = sum(y(:, idx).');
  else
   modes(:, i) = y(:, idx).';
  end
 end
end