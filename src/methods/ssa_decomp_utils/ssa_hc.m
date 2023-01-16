function modes = ssa_hc(x, L, nc, epsilon)
% modes = ssa_hc(x, L, nc, epsilon)
%
% Apply SSA combined with Hierarchical Clustering for automatic components
% extraction and grouping
%
%  INPUT:
%  x:       input time series
%  L:       Window length (embedding dimension)
%  nc:      Number of components to extract
%  epsilon: Component energy contribution threshold (default 3e-2) the last components whose energy contribution is relatively
% null are eliminated
%
%  OUTPUT:
%  modes:   extracted modes
%
% Ref: [J. Harmouche, D.Fourer, F. Auger, P.Flandrin and P. Borgnat. The Sliding Singular Spectrum Analysis:
% a Data-Driven Non-Stationary Signal Decomposition Tool. IEEE TSP 2016.]

if ~exist('epsilon', 'var')
 epsilon = 3e-2; %0.5
end

Nx=length(x);
 
%modes = zeros(Nx, 1);
if nc < 1
 %modes = zeros(Nx, 1);
 return;
end

[first_y, d] = ssa(x,L,1);


% plot(d)
% 
% pause

% get the number of non negligible singular values to consider
supk = length(find( (d/max(d))>epsilon));

if supk < nc   %% the number of singular values should be at least equal to nc
 warning('epsilon is too high, trying with eps=%.2f', d(nc)/max(d))
 supk = nc;
end

%% extract the other elementary components
y = zeros(Nx, supk);
y(:,1) = first_y;
for I = 2:supk;
 [y(:,I),~]=ssa(x,L,I); 
end
   
%wcorr = pdist(x.', 'Euclidean');
wcorr = 1-squareform(pdist(y.', 'correlation'));
%wcorr = wCorrMat2(x,L,supk);

%pause

if ~isempty(wcorr)
 Z = linkage(wcorr);                       
 T = cluster(Z,'maxclust', nc);  
 
 %dendrogram(Z)
 %title(sprintf('Number of components to extract: %d', nc));
 %pause
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
