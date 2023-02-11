function [xr,mask_total] = pb_method(x, Ncomp, use_sst, ds, beta, alpha, div, Pnei, PneiMask, M, L,return_comps, return_instf)
% This is the function that is called from the Python-based benchmark.


% Get signal parameters and compute the TF representation.
N = length(x);

% Bayesian method default parameters
if nargin<2 || isempty(Ncomp)
    Ncomp = 1;
end

if nargin<3 || isempty(use_sst)
    use_sst = false;
end

if nargin<4 || isempty(ds)
    ds    = 3;    % variance of the random walk in the temporal model
end

if nargin<5 || isempty(beta)
    beta  = 0.6;
end

if nargin<6 || isempty(alpha)
    alpha = 0.5;
end

if nargin<7 || isempty(div)
    div   = 4;                         % 1 = KL
end

if nargin<8 || isempty(Pnei)
    Pnei = 15;
end

if nargin<9 || isempty(PneiMask)
    PneiMask = 15;
end

if nargin<10 || isempty(M)
    M = N;
end

if nargin<11 || isempty(L)
    L = 10;
end

if nargin<12 || isempty(return_comps)
    return_comps = false;
end

if nargin<13 || isempty(return_instf)
    return_instf = false;
end

detect = 0;
ifplot =  0;

% Tfr parameter
%L = 30;
%M = 500;
[tfr,stfr]  = tfrsgab2(x, M, L);
reconstruct_func = @(atfr,L,M)rectfrgab(atfr, L, M);

% Use synchrosqueezed version
if use_sst
    tfr = stfr;
    reconstruct_func = @(atfr,L,M)rectfrsgab(atfr, L, M);
end

% Apply the method and generate a mask per component
[mask,instf] = pseudoBay(tfr, Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect, PneiMask);

% Generate a combined mask of all components.
mask_total = sum(mask,3);
mask_total(mask_total~=0) = 1;

% Recover components and signal
x_hat = zeros(N,Ncomp);

% Inversion of the masked STFT.
for c = 1:Ncomp
    x_hat(:,c) = real(reconstruct_func(tfr .* mask(:,:,c), L, M));
    %     [xr,~] = tfristft(tfr .* mask(:,:,c),1:N,w,0);
end

% Return reconstructed signal.
% xr = real(reconstruct_func(tfr .* mask_total, L, M));
xr = sum(x_hat,2).';

if return_comps
    xr = x_hat.';
end

if return_instf
    xr = instf.'/M;
end