clear                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          % clear all
close all
% 
% 
%folder = './';
%addpath(strcat([folder 'tools']));
% addpath(strcat([folder 'mfiles']));
% addpath(strcat([folder 'synchrosqueezedSTFT']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length
% X(:,1)    = (fmconst(N, 0.1));
% X(:,1) = fmlin(N,0.2,0.4);


%X(:,1) = fmsin(N,0.1,0.3,320,1,0.3,+1);


X(:,1) = 0.8 * fmlin(N,0.41,0.1);
X(:,2) = fmlin(N,0.1,0.45);

x0 = sum(X, 2);

x0 = sigmerge(x0, randn(size(x0)), 10);

%% Ground truth
[tfr]  = tfrgab2(x0, M, L);
Y0 = abs(tfr(1:M/2,:)).^2;
imagesc(abs(tfr))

%%
nb_comp = 2;

tmp = abs(tfr(1:M/2,:)).^2;

mask  = zeros(size(tmp,1),size(tmp,2), nb_comp);
y_hat = zeros(size(X));

seuil = 1e-2;

for i = 1:nb_comp
    
 Y = tmp;
 Ns = 1;
 [Fct]=comp_Fc(M,L);                  %% Data distribution
 Fct = repmat(Fct,[1,1,Ns]);

%% Initialize Parameters
c=1e-4; % TV regularization parameter
cl=1e-2; % Laplacian spatial prior reg hyperparameter
step_Nx = 1; % Depth grid subsampling
stepg = 1e-3; % Gradient step for the M-step
ifplot = 0;  %%debug

%% Prior choice
% reg_mu='None';
% reg_mu='Lap';
reg_mu='TV';

%% EM algorithm
tic
% [a_out,T_out]=EM_ranging(Y',Fc,Ncomp,c,si,step_Nx,stepg,nb_label);
[W_out,P,ind0]=Mod_Estim_W_EM(Y',Fct,Ns,step_Nx,reg_mu,c,cl,ifplot);

toc
%% Plot reconstructions

t_axis = (1:N)-1;
Mh = round(M/2);
f_axis = ((1:Mh)-1)/M;

xx = zeros(size(tmp));
xx(P.'> seuil) = 1;
mask(:,:,i) = xx;

tmp(P.'> seuil) = 0;   %% remove component

end



