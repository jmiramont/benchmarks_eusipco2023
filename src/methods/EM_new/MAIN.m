clear                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     % clear all
close all
clc
warning('off','all')

folder = './';
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'synchrosqueezedSTFT']));


%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length
% X(:,1)    = (fmconst(N, 0.1));
% X(:,1) = fmlin(N,0.2,0.4);
X(:,1) = fmsin(N,0.3,0.45,320,1,0.3,+1);

x0 = sum(X,2);
Ncomp = size(X,2);                  %% number of components

SNR = -5;
% add noise
x = sigmerge(x0, randn(size(x0)), SNR);

%% Ground truth
[tfr]  = tfrgab2(x0, M, L);
Y0 = abs(tfr(1:M/2,:)).^2;

%% Data
[tfr]  = tfrgab2(x, M, L);
Y = abs(tfr(1:M/2,:)).^2;


[Fc]=comp_Fc(M,L);                  %% Data distribution
Fc = repmat(Fc,[1,1,Ncomp]);

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
[W_out,P,ind0]=Mod_Estim_W_EM(Y',Fc,Ncomp,step_Nx,reg_mu,c,cl,ifplot);
toc
%% Plot reconstructions

t_axis = (1:N)-1;
Mh = round(M/2);
f_axis = ((1:Mh)-1)/M;

figure(111)
imagesc(t_axis, f_axis,  Y);
set(gca, 'Ydir', 'normal');
xlabel('time [samples]')
ylabel('normalized frequency')

figure(222)
imagesc(t_axis, f_axis,  P.');
set(gca, 'Ydir', 'normal');
xlabel('time [samples]')
ylabel('normalized frequency')
%figure(5)
% hold on
% plot(amp,'b')
% plot(W_out,'r')
% hold off
% legend('GT','est')
