clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative robustness evaluation for mask estimation
%  associated to signals whose frequency component are not defined for all
%  time.
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 23-mar-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'SSA']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Compute_Amplitude_DF']));


snr_range = [-20 20]; % SNR range to compute
MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 1;                          %% number of components
N     = 500;                        %% signal length
MD0 = ones(N,1);
MD0(1:200,1)   = 0;   % Comp 1
MD0(401:500,1) = 0;   % Comp 1
X(:,1) = MD0.*(fmconst(N, 0.1));
x0 = sum(X,2);
X = transpose(X);
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 10; % variance of the random walk in the temporal model
n_pad = 15; 
detect=1;
div = 4;

alpha = 0.5;
beta = 0.5;
%% Compute ground truth
tf0=zeros(N,Ncomp);
Lbest=20;
[S]  = tfrgab2(fmconst(N, 0.1), M, Lbest);
S=(abs(S(1:round(M/2),:)));
Nexist = find(MD0~=0);
lN = length(Nexist);
for i=1:lN
    ii = Nexist(i);
    [~,mm]=max(S(:,ii));
    tf0(ii)=mm;
end
tf0 = max(tf0,1);
[tfr]  = tfrgab2(x0, M, L);
Sp = abs(tfr(1:M/2,:)).^2;



%% Initialization
methods_name = {'Proposed-\alpha=0.4,\beta=0.4,sAB-div.',...
                'Proposed-\alpha=0.2,\beta=0.4,sAB-div.',...
                'Proposed-\alpha=0.4,\beta=0.2,sAB-div.',...
                'Proposed-\alpha=0.2,\beta=1.2,sAB-div.',...
                'Proposed-\alpha=0.7,\beta=1.2,sAB-div.',...
                'Proposed-\alpha=0.2,\beta=1.5,sAB-div.',...
                'Oracle'
                };
            
            
methods_to_use = [1 2 3 4 5 6 7];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);

figure
colormap(flipud(gray))
subplot(1,2,1)
imagesc(Sp)
axis square
ax = gca;
ax.YDir = 'normal'
xlabel('Time index','FontSize', 10, 'FontWeight', 'bold')
ylabel('Normalized frequency','FontSize', 10, 'FontWeight', 'bold')
yticks([1 M/4 M/2])
yticklabels({0,0.25,0.5})

subplot(1,2,2)
plot(MD0)
axis square
%ax = gca;
%ax.YDir = 'normal'
xlabel('Time index','FontSize', 10, 'FontWeight', 'bold')
ylabel('Presence array','FontSize', 10, 'FontWeight', 'bold')
ylim([-0.25 1.25])
xlim([1 N])
%yticklabels({-0.25,1.25})