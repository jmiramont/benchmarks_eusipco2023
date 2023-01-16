function signal_r = method_ssa_decomp(signal, n_components, L, epsilon, return_comps)
%METHOD_SSA_DECOMP Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3 || isempty(L)
L  = 40;  %% embedded dimension SSA parameter
end

if nargin < 4 || isempty(epsilon)
epsilon = 5e-3;    %% singular spectrum thresholding parameter (increase for more robustness to noise)
end

if nargin < 5 || isempty(return_comps)
    return_comps = false;
end

Y = ssa_decomp(signal, L, n_components, epsilon);

signal_r = sum(Y,2).';

if return_comps
    signal_r = Y.';
end

end

