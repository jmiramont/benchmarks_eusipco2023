function xr = em_method(x,Ncomp,M,L,c,cl,step_Nx,stepg,seuil,return_comps)
% This function is called from the (python based) benchmark.
% It wraps the EM method and parameters (this of course can be modified).
%--------------------------------------------------------------------------

if ~iscolumn(x)
    x = x.';
end
N = length(x);

%% Define the default parameters -------------------------------------------
if nargin<3 || isempty(M)
    M       = N;       %% nombre de bins frequentiels
end

if nargin<4 || isempty(L)
    L       = 20;        %% taille de la fenetre d'analyse en samples
end

if nargin<5 || isempty(c)
    c=1e-4; % TV regularization parameter
end

if nargin<6 || isempty(cl)
    cl=1e-2; % Laplacian spatial prior reg hyperparameter
end

if nargin<7 || isempty(step_Nx)
    step_Nx = 1; % Depth grid subsampling
end

if nargin<8 || isempty(stepg)
    stepg = 1e-4; % Gradient step for the M-step
end

if nargin<9 || isempty(seuil)
    seuil = 1e-7;
end

% Return components.
if nargin<10 || isempty(return_comps)
    return_comps = false;
end

%% Apply EM method and get the masks --------------------------------------
x_hat = zeros(N,Ncomp);
[tfr,stfr]  = tfrsgab2(x, M, L);
Y = abs(tfr(1:M/2,:)).^2;
% [ y_hat, mask ] = EM_estim_multi(tfr, L, Ncomp, seuil, c, cl, step_Nx, 0);

Ns = 1;
[Fct]=comp_Fc(M,L);                  %% Data distribution
Fct = repmat(Fct,[1,1,Ns]);
ifplot = 0;  %%debug

%% Prior choice
% reg_mu='None';
% reg_mu='Lap';
reg_mu='TV';

%% EM algorithm
for c = 1:Ncomp
    % [a_out,T_out]=EM_ranging(Y',Fc,Ncomp,c,si,step_Nx,stepg,nb_label);
    [W_out,P,ind0]=Mod_Estim_W_EM(Y',Fct,Ns,step_Nx,reg_mu,c,cl,ifplot);

    xx = zeros(size(Y));
    xx(P.'> seuil) = 1;
    mask = [xx;xx(end:-1:1,:)];
    
    Y(P.'> seuil) = 0;   %% remove component

    x_hat(:,c) = real(rectfrgab(tfr .* mask, L, M));
end

% Generate a combined mask of all components and invert the masked STFT.
% mask_total = sum(mask,3);
% mask_total(mask_total~=0) = 1;
% mask_total = ones(size(tfr))
% xr = real(rectfrgab(tfr .* mask_total, L, M));


xr = sum(x_hat,2);

% If return_comps is True, then return the components insted of the signal.
if return_comps
    xr = x_hat.';
end
