clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the comparative robustness evaluation for mask estimation only
%  using the RMSE to reproduce the results in Fig 2 of the journal paper
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
Ncomp = 1;                          %% number of components
N     = 500;                        %% signal length
x0    = real(fmlin(N,0.13,0.3));    %% linear chirp
n_pad = 15;
%% Method parameters
Pnei = 0; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 0;
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 2; % variance of the random walk in the temporal model
alpha = 0.5;
beta = 0.5;
%% Compute ground truth
tf0=zeros(N,1);

[tfr]  = tfrgab2(x0, M, L);
tfr=(abs(tfr(1:round(M/2),:)));
for i=1:N
    [~,mm]=max(tfr(:,i));
    tf0(i)=mm;
end
detect=0;





 

%% Initialization
% methods_name = {'Brevdo-SST',                      'Brevdo-STFT',...
%                 'Proposed-\beta=0.3,Beta-div.',    'Proposed-\beta=0.7,Beta-div.',...
%                 'Proposed-\beta=1,Beta-div.',      'Proposed-KL-div.',...
%                 'Proposed-\alpha=0.3,Renyi-div.',  'Proposed-\alpha=0.5,Renyi-div.',...
%                 'Proposed-\alpha=0.8,Renyi-div.'      
%                 };
methods_name = {'PB-\alpha=0.5',...
                'PB-\alpha=0.8',...
                'PB-\beta=0.7',...
                'PB-\beta=0.9',...
                'Proposed-\alpha=0.4,\beta=0.4,sAB-div.',...      % Robust            - Mass covering  -   5 - old
                'Proposed-\alpha=0.2,\beta=1.5,sAB-div.',...      % Outlier focussing - Mass covering  -   6 - old
                'Proposed-\alpha=1.2,\beta=-0.4,sAB-div.',...     % Robust            - Mode seeking   -   7 - old
                'Proposed-\alpha=1.1,\beta=1.1,sAB-div.',...      % Outlier focussing - Mode seeking   -   8 - old
                };
            
methods_to_use = [1 2 3 4 5 6 7 8];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);

L2ErPos_out = zeros(length(SNRt), nb_methods);


%% Compute RMSE
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        L2ErPos_tmp = zeros(1,length(MCrep));

        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            

            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            
            switch(methods_to_use(ind_met))
                case 1  %% Rényi divergence
                        alpha  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 3;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,PneiMask);

                case 2  %% Rényi divergence
                        alpha  = 0.8; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 3;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,PneiMask);

                case 3  %% Beta divergence
                        beta  = 0.7; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 2;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,PneiMask);

                case 4  %% Beta divergence
                        beta  = 0.9; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 2;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);

                case 5 %% Alpha-Beta divergence
                        beta  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 4;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);

                case 6  %% Alpha-Beta divergence
                        beta  = 1.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha  = 0.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 4;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);

                case 7 %% Alpha-Beta divergence
                        beta  = -0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha  = 1.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 4;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);

                case 8  %% Alpha-Beta divergence
                        beta  = 1.1; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        alpha  = 1.1; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 4;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,PneiMask);

            end  %% switch

            % Updated MSE - not relative anymore
            L2ErPos_tmp(it) = 10*log10(sum((tf(n_pad:(end-n_pad))-tf0(n_pad:(end-n_pad))).^2)./(N*M*M));
        end    %% for methods
        L2ErPos_out(indsnr, ind_met) = mean(L2ErPos_tmp);
    end  %% methods
 
end %% snrs
   

%% Plot
% cols         = {'k-x'  'b-o' 'g-s' 'b-.' 'g-.' 'r-.' 'r-v'   'b-.' 'b--' 'b*'};
cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-o' 'g-o' 'r-o' 'k--' 'b--' 'g--' 'r--' 'k-v' 'b-v' 'g-v'};
% cols         = {'k-x'  'b-x' 'g-x' 'r-x' 'k-o'  'b-o' 'g-o' 'r-o' 'k-v'  'b-v' 'g-v' 'r-v'};
leg = {};

figure(1)
for ind_met =  1:nb_methods
    
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, L2ErPos_out(:,ind_met), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 14)
 ylabel('RMSE (dB)', 'FontSize', 14)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

grid


