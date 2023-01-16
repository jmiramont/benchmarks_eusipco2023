clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compare the estimation performance of competings methods when using
%  different non binary mask to perform signal synthesis for varying SNR.
%  This code generates figures 8, 9 and 10 of the journal paper.
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 31-may-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));



snr_range = [-20 20]; % SNR range to compute
MCrep = 50;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 3;                          %% number of components
N     = 500;                        %% signal length
X(:,1) = real(fmconst(N, 0.1));
X(:,2) = real(fmlin(N,0.13,0.3));
X(:,3) = real(fmsin(N,0.3,0.45,320,1,0.3,+1));

x0 = sum(X,2);
X = transpose(X);

[tfr]  = tfrgab2(x0, M, L);
Spect = abs(tfr(1:M/2,:)).^2;

figure
colormap(flipud(gray))
imagesc(Spect)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 12, 'FontWeight', 'bold')
ylabel('Normalized frequency','FontSize', 12, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})