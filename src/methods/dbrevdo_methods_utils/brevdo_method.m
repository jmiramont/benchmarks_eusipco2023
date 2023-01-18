function xr = brevdo_method(x, Ncomp, use_sst, Pnei, M, L, return_comps, return_freq)
%% This is the function that is called from the Python-based benchmark.
N = length(x);

% Default parameters
if nargin <3 || isempty(use_sst)
    use_sst=true;
end

if nargin<4 || isempty(Pnei)
    Pnei=35;
end

if nargin<5 || isempty(M)
    M=N;
end

if nargin<6 || isempty(L)
    L=10;
end

if nargin<7 || isempty(return_comps)
    return_comps = false;
end

if nargin<8 || isempty(return_freq)
    return_freq = false;
end

% Time-frequency representation parameters
% M       = 500;       %% Number of frequential bin
% L       = 20;        %% analysis window size (in bin)

% % Tfr parameter
% L = 30;
% M = 500;
[tfr,stfr]  = tfrsgab2(x, M, L);
reconstruct_func = @(atfr,L,M)rectfrgab(atfr, L, M);

% Use synchrosqueezed version
if use_sst
    tfr = stfr;
    reconstruct_func = @(atfr,L,M)rectfrsgab(atfr, L, M);
end

% Apply the method
[~, mask, Cs] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);

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
xr = real(reconstruct_func(tfr .* mask_total, L, M));
% xr = sum(x_hat,2).';

% Return the individual components.
if return_comps
    xr = x_hat.';
end

% Return the individual inst freq.
if return_freq
    xr = Cs/M;
end

