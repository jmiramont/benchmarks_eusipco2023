clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.1 of the paper ESTIMATION OF OVERLAPPING MODES USING SPARSE MODELING OF SIGNALINNOVATION
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 16-feb-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'PseudoBay']));


snr_range = [-20 20]; % SNR range to compute
MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 3;                          %% number of components
N     = 500;                        %% signal length

%%%%%%%%%%%%%%%%%%%%%%% non overlaping signal
X1 = zeros(N,Ncomp);
X1(:,1) = (fmconst(N, 0.1));
X1(:,2) = (fmlin(N,0.13,0.3));
X1(:,3) = (fmsin(N,0.3,0.45,320,1,0.3,+1));

x10 = sum(X1,2);
X1 = transpose(X1);

[tfr1] = tfrgab2(x10, M, L);       %% compute STFT
Spect1 = abs(tfr1(1:M/2,:)).^2;

%%%%%%%%%%%%%%%%%%%%%%% overlaping signal
X2 = zeros(N,Ncomp);
X2(:,1) = (fmlin(N,0.05,0.45));
X2(:,2) = (fmlin(N,0.4,0.15));
X2(:,3) = (fmsin(N,0.15,0.4,320,1,0.15,+1));

x20 = sum(X2,2);
X2 = transpose(X2);

[tfr2] = tfrgab2(x20, M, L);       %% compute STFT
Spect2 = abs(tfr2(1:M/2,:)).^2;


%%%%%%%%%%%%%%%%%%%%%%% Plot


figure
colormap(flipud(gray))
subplot(1,2,1)
imagesc(Spect1)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 10, 'FontWeight', 'bold')
ylabel('Normalized frequency','FontSize', 10, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})

subplot(1,2,2)
imagesc(Spect2)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 10, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})








