function [xr,tf2,mask_total] = em_method(x,Ncomp,M,L, Pnei, step_r, step_v, return_comps, return_freq, bolpol)
% This function is called from the (python based) benchmark.
% It wraps the EM method and parameters (this of course can be modified).
%--------------------------------------------------------------------------

if ~iscolumn(x)
    x = x.';
end
N = length(x);

%% Define the default parameters -------------------------------------------
if ~exist('M', 'var') || isempty(M)
    M       = N;       %% nombre de bins frequentiels
end

if ~exist('L', 'var') || isempty(L)
    L       = 20;        %% taille de la fenetre d'analyse en samples
end

if ~exist('c', 'var') || isempty(c)
    c=1e-4; % TV regularization parameter
end

if ~exist('Pnei', 'var') || isempty(Pnei)
    Pnei=48;
end

if ~exist('step_r', 'var') || isempty(step_r)
    step_r= 3*ceil(sqrt(M/(pi*L)))+1; 
end

if ~exist('step_v', 'var') || isempty(step_v)
    step_v= 4*step_r;
end

% if ~exist('cl', 'var')
%     cl=1e-2; % Laplacian spatial prior reg hyperparameter
% end
% if ~exist('step_Nx', 'var')
%     step_Nx = 1; % Depth grid subsampling
% end
% 
% if ~exist('seuil', 'var')
%     seuil = 1e-7;
% end

% Return components.
if ~exist('return_comps', 'var') || isempty(return_comps)
    return_comps = false;
end

% Return components.
if ~exist('return_freq', 'var') || isempty(return_comps)
    return_freq = false;
end

if ~exist('bolpol', 'var') || isempty(bolpol)
    bolpol=1;
end

ifplot = 0;
% Parameter for sequential MMAP estimation
% step_r = 30; %30; % removal window size
% step_v = 1;  %1 neighborhood search size

%% Prior choice
% reg_mu='None';
reg_mu='Lap';
%reg_mu='TV';

%%zero-padding
zp = 0; %round(3*L); %%3 sigma thumb rule 
z = zeros(zp,1);
z0 = zp+1;
z1 = z0+N-1;

%% Apply EM method and get the masks --------------------------------------
x_hat = zeros(N+2*zp,Ncomp);

M2 = floor(M/2);
[tfr]  = tfrgab2([z;x;z], M, L);
Spect = abs(tfr(1:M2,:)).^2;


% Ns = 1;
[Fct]=(comp_Fc(M,L))+eps; %% Data distribution
[W_out,P,tf]=Mod_Estim_W_EM_multi(Spect',Fct,Ncomp,1,reg_mu,c,step_r,step_v,ifplot,bolpol);


%% EM algorithm mask
delta_m = Pnei;
mask_total = 0;
for c = 1:Ncomp

    xx = zeros(size(Spect));
    
    for n = 1:N
        m0 = max(1,tf(n,c)-delta_m);
        m1 = min(tf(n,c)+delta_m, M2);
        
        xx(m0:m1,n) = c; %W_out(n);
    end
    
    mask = [xx;xx(end:-1:1,:)];
    mask_total = mask_total + mask; % Combine all masks to compare later.
    %% debug
%      figure;
%      imagesc(abs(tfr) .* mask)
%      pause
    x_hat(:,c) = real(rectfrgab(tfr .* logical(mask), L, M));
end
x_hat = x_hat(z0:z1,:);

% Generate a combined mask of all components and invert the masked STFT to
% obtain a reconstruction of the signal.
% mask_total(mask_total~=0) = 1;
xr = real(rectfrgab(tfr .* logical(mask_total), L, M));

% Sum the components to generate the output signal -> JMM 27/02: This is a
% problem when the individual masks overlap!
% xr = sum(x_hat,2).';

% If return_comps is True, then return the components insted of the signal.
if return_comps
    xr = x_hat.';
end

% If return_freq is True, then return the components' IF
tf2 = tf.'/M;
if return_freq
    xr = tf.'/M;
end