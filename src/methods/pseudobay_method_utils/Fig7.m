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
PneiMask = 10;
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

MAE_out = zeros(length(SNRt), nb_methods);

%% Compute RQF
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        MAE_tmp = zeros(Ncomp,MCrep);
        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
            
            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            switch(methods_to_use(ind_met))
                case 1  %% Beta divergence
                        %% Bayesian method parameters
                        alpha  = 0.4;
                        beta   = 0.4;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                case 2  %% Renyi diverg
                        %% Bayesian method parameters
                        alpha  = 0.2;
                        beta   = 0.4;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                case 3  %% Beta divergence
                        alpha  = 0.4;
                        beta   = 0.2;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                case 4  %% Beta divergence
                        %% Bayesian method parameters
                        alpha  = 0.2;
                        beta   = 1.2;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                case 5  %% Beta divergence
                        %% Bayesian method parameters
                        alpha  = 0.7;
                        beta   = 1.2;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                case 6  %% Renyi diverg
                        %% Bayesian method parameters
                        alpha  = 0.2;
                        beta   = 1.5;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,~,MD] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                 case 7 %% Oracle
                        %% Bayesian method parameters
                        [tfr]  = tfrgab2(x, M, L);
                        [MD] = detectOracle(tfr,M,L,tf0);
                        
            end  %% switch

            % MAE
            for Nc = 1:Ncomp
                MAE_tmp(it) = sum(abs(MD - MD0));
            end
            
            
%             figure
%             subplot(2,1,1)
%             imagesc(abs(tfr.^2))
%             subplot(2,1,2)
%             imagesc(sum(mask(1:M/2,:,:),3))
                        
        end    %% for methods
        MAE_out(indsnr, ind_met) = mean(MAE_tmp);
    end  %% methods
 
end %% snrs
   
for Nc = 1:Ncomp
    %% Plot
    cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-s' 'g-.' 'r-.' 'b-v'  'b--' 'b*'};
    leg = {};

    figure(Nc)
    for ind_met =  1:nb_methods

     if ind_met == 1
      hold off
     else
      hold on
     end

     h(ind_met) = plot(SNRt, squeeze(MAE_out(:,ind_met)), cols{ind_met});
     xlabel('SNR (dB)', 'FontSize', 14)
     ylabel('MAE', 'FontSize', 14)
     leg{ind_met} = methods_name{methods_to_use(ind_met)};
    end
    legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

    grid
end

