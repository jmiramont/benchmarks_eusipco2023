clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Comparative study of ridge extraction in the presence of overlapping
%  ridges
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 16-feb-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'SSA']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'PseudoBay']));


snr_range = [-20 20]; % SNR range to compute
MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)
%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 1;
N=500;
X = zeros(N,Ncomp);
X(:,1) = real(fmconst(N, 0.1));

x0 = sum(X,2);
X = transpose(X);
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 20;
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 1; % variance of the random walk in the temporal model
detect = 0;
n_pad = 15;

%% Bayesian method parameters
beta  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
alpha = 0.5; % Renyi divergence hyperparameter ||  POSITIVE AND DIFFERENT TO 1
div   = 3;   % 1 = KL
             % 2 = beta
             % 3 = Renyi

             
sigmt = [0.1,1,2.5,5,10,100];
% sigmt = [50,100];
SNRt = snr_range(1):2:snr_range(2);
RQF_out = zeros(Ncomp,length(SNRt),length(sigmt));

%% Compute RQF
for indsnr = 1:length(SNRt)
	fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
	SNRi = SNRt(indsnr);
    
    for indsigm = 1:length(sigmt)
        sigm_conv = sigmt(indsigm);
    
        RQF_tmp = zeros(Ncomp,MCrep);
        for it = 1:MCrep   %% iterations
                clc;
                disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
                disp(strcat(['Sigma :', num2str(sigm_conv)]));
                disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
                % add noise
                x = sigmerge(x0, randn(size(x0)), SNRi);

                [tfr]  = tfrgab2(x, M, L); % compute tfr

                
                [mask,~,~,~,~,~] = pseudoBay_conv(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,sigm_conv,PneiMask);
%                 [mask,~,~,~,~,~] = pseudoBay_filt(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect);
                
                x_hat = zeros(Ncomp,N);
                for c = 1:Ncomp
                   x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));       
                end
                % Match components and reordering for comparison
%                 [I,~] = match_components(X, x_hat); 
%                 x_hat = x_hat(I,:);

%                 figure(10)
%                 hold on
%                 plot(real(X(n_pad:end-n_pad)),'b')
%                 plot(real(x_hat(n_pad:end-n_pad)),'r')
%                 hold off

                for Nc = 1:Ncomp
                    RQF_tmp(Nc,it) = RQF(real(X(Nc,n_pad:end-n_pad)), x_hat(Nc,n_pad:end-n_pad));
                end
%                 RQF_tmp(Nc,it)
        end
        RQF_out(:,indsnr,indsigm) = mean(RQF_tmp,2);
    end
end


for Nc = 1%:Ncomp
    %% Plot
    cols         = {'r--' 'g-o' 'r-s' 'k--' 'b-*' 'k-s'};
    leg = {};

    figure(Nc)
    for indsigm = 1:length(sigmt)

     if indsigm == 1
      hold off
     else
      hold on
     end

     h(indsigm) = plot(SNRt, squeeze(RQF_out(Nc,:,indsigm)), cols{indsigm});
     xlabel('SNR (dB)', 'FontSize', 14)
     ylabel('RQF (dB)', 'FontSize', 14)
     leg{indsigm} = strcat(['Std :', num2str(sigmt(indsigm))]);
    end
    legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

    grid
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


