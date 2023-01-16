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
MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 3;                          %% number of components
N     = 500;                        %% signal length
X(:,1) = real(fmconst(N, 0.1));
X(:,2) = real(fmlin(N,0.13,0.3));
X(:,3) = real(fmsin(N,0.3,0.45,320,1,0.3,+1));

%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 20;
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 1; % variance of the random walk in the temporal model
n_pad = 15; 

%% Compute ground truth
% X(1:200,1)   = 0;   % Comp 1
% X(401:500,1) = 0;   % Comp 1
% X(251:300,2) = 0;   % Comp 2

detect=0;


x0 = sum(X,2);
X = transpose(X);

[tfr]  = tfrgab2(x0, M, L);
Spect = abs(tfr(1:M/2,:)).^2;

figure
imagesc(Spect)
axis square
xlabel('Time')
ylabel('Normalized frequency')
yticks([1 125 250])
yticklabels({'0','0.25','0.5'})

%% Initialization
methods_name = {'Brevdo-STFT',...
                'Proposed-\alpha=0.4,\beta=0.4,sAB-div.',...
                'SST-\alpha=0.4,\beta=0.4,sAB-div',...
                'NB-\sigma=5-\alpha=0.4,\beta=0.4,sAB-div.',...
                'NB-\sigma=7-\alpha=0.4,\beta=0.4,sAB-div.',...
                'NB-\sigma=10-\alpha=0.4,\beta=0.4,sAB-div.',...
                'NB-\sigma=20-\alpha=0.4,\beta=0.4,sAB-div.'
                };
            
beta  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
alpha  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0

methods_to_use = [1 2 3 4 5 6 7];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);


RQF_out = zeros(Ncomp,length(SNRt), nb_methods);





%% Compute RQF
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        RQF_tmp = zeros(Ncomp,MCrep);
        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
            
            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            
            switch(methods_to_use(ind_met))
                case 1  %% Brevdo STFT
                        [tfr] = tfrgab2(x, M, L);       %% compute STFT
                        [~, mask] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                           x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));       
                        end
                case 2  %% Alpha-Beta divergence
                        div   = 4;   % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                            x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
                        end
                 case 3 %% SST Alpha-Beta divergence
                        div   = 4;    % 3 = Renyi
                        [tfr,stfr] = tfrsgab2(x, M, L); %% compute SST
                        [mask,~,~,~,~,~] = pseudoBay_sst(stfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                            x_hat(c,:) = real(rectfrsgab(stfr .* mask(:,:,c), L, M));
                        end
                case 4 %% Non-binary mask Alpha-Beta divergence
                        div   = 4;    % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        sigm_conv = 5;
                        [mask,~,~,~,~,~] = pseudoBay_conv(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,sigm_conv,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                            x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
                        end
                case 5 %% Non-binary mask Alpha-Beta divergence
                        div   = 4;    % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        sigm_conv = 7;
                        [mask,~,~,~,~,~] = pseudoBay_conv(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,sigm_conv,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                            x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
                        end
                case 6  %% Non-binary mask Alpha-Beta divergence
                        div   = 4;    % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        sigm_conv = 10;
                        [mask,~,~,~,~,~] = pseudoBay_conv(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,sigm_conv,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                             x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
                        end
                case 7  %% Non-binary mask Alpha-Beta divergence
                        div   = 4;    % 3 = Renyi
                        [tfr]  = tfrgab2(x, M, L);
                        sigm_conv = 20;
                        [mask,~,~,~,~,~] = pseudoBay_conv(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot,detect,sigm_conv,PneiMask);
                        x_hat = zeros(Ncomp,N);
                        for c = 1:Ncomp
                             x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
                        end

                        
            end  %% switch

            % Match components and reordering for comparison
            [I,~] = match_components(X, x_hat); 
            x_hat = x_hat(I,:);

            % RQF
            for Nc = 1:Ncomp
                RQF_tmp(Nc,it) = RQF(real(X(Nc,n_pad:end-n_pad)), x_hat(Nc,n_pad:end-n_pad));
            end
            
        end    %% for methods
        RQF_out(:,indsnr, ind_met) = mean(RQF_tmp,2);
    end  %% methods
 
end %% snrs
   

for Nc = 1:Ncomp
    %% Plot
    cols         = {'r-o' 'k-o' 'g-o' 'b-x' 'r-s' 'k-s' 'b-x' 'r-.' 'r-v'   'r-o'  'k--' 'k-*'};

    
    leg = {};

    figure(Nc)
    for ind_met =  1:nb_methods

     if ind_met == 1
      hold off
     else
      hold on
     end
     
     h(ind_met) = plot(SNRt, squeeze(RQF_out(Nc,:,ind_met)),cols{ind_met});
     xlabel('SNR (dB)', 'FontSize', 14)
     ylabel('RQF (dB)', 'FontSize', 14)
     leg{ind_met} = methods_name{methods_to_use(ind_met)};
    end
    legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

    grid
end
