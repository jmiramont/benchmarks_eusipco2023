clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compare the reconstruction performance of the competing methods in term 
%  of the RQF with a SNR = 10dB to reproduce the results depicted in Table 1
%  of the journal paper
%  
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 31-may-2021



folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));



snr_range = [-20 20]; % SNR range to compute
MCrep = 20;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

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
n_pad = 15;
%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 10;
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 2; % variance of the random walk in the temporal model
alpha = 0.5;
beta = 0.5;
detect = 0;
%% Compute ground truth
tf0=zeros(N,1);

[tfr]  = tfrgab2(x0, M, L);
tfr=(abs(tfr(1:round(M/2),:)));
for i=1:N
    [~,mm]=max(tfr(:,i));
    tf0(i)=mm;
end

 

%% Initialization
% methods_name = {'Brevdo-SST',                      'Brevdo-STFT',...
%                 'Proposed-\beta=0.3,Beta-div.',    'Proposed-\beta=0.7,Beta-div.',...
%                 'Proposed-\beta=1,Beta-div.',      'Proposed-KL-div.',...
%                 'Proposed-\alpha=0.3,Renyi-div.',  'Proposed-\alpha=0.5,Renyi-div.',...
%                 'Proposed-\alpha=0.8,Renyi-div.'      
%                 };
methods_name = {'Brevdo-STFT',...
                'Proposed-\alpha=0.4,\beta=0.4,sAB-div.',...      % Outlier focussing - Mode
                'Proposed-\alpha=0.2,\beta=0.4,sAB-div.',...      % Outlier focussing - Mode
                'Proposed-\alpha=0.4,\beta=0.2,sAB-div.',...      % Outlier focussing - Mode
                'Proposed-\alpha=0.2,\beta=1.2,sAB-div.',...      % Robust            - Mode 
                'Proposed-\alpha=0.7,\beta=1.2,sAB-div.',...      % Robust            - Mode 
                'Proposed-\alpha=0.2,\beta=1.5,sAB-div.',...      
                };
            
methods_to_use = [1 2 3 4 5 6 7];   % insert here the indices of the methods to compare (names above)
nb_methods = length(methods_to_use);


SNR = 10;

RQF_out = zeros(Ncomp, nb_methods);
RQF_var_out = zeros(Ncomp, nb_methods);

%% Compute RMSE
for ind_met = 1:length(methods_to_use)
    L2ErPos_tmp = zeros(1,length(MCrep));


    for it = 1:MCrep   %% iterations
        clc;
        disp(strcat(['SNR : ',num2str(SNR)]));
        disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
        disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))


        % Add noise
        x = sigmerge(x0, randn(size(x0)), SNR);


        switch(methods_to_use(ind_met))
             case 1  %% Brevdo STFT
                    [tfr] = tfrgab2(x, M, L);       %% compute STFT
                    [~, mask] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                       x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));       
                    end
             case 2 %% Alpha-Beta divergence
                    beta  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M)); 
                    end
            case 3  %% Alpha-Beta divergence
                    beta  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));      
                    end
            case 4 %% Alpha-Beta divergence
                    beta  = 0.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.4; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M));      
                    end
            case 5  %% Alpha-Beta divergence
                    beta  = 1.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M)); 
                    end
            case 6 %% Alpha-Beta divergence
                    beta  = 1.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.7; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
                    x_hat = zeros(Ncomp,N);
                    for c = 1:Ncomp
                        x_hat(c,:) = real(rectfrgab(tfr .* mask(:,:,c), L, M)); 
                    end
            case 7  %% Alpha-Beta divergence
                    beta  = 1.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    alpha  = 0.2; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                    div   = 4;   % 3 = Renyi
                    [tfr]  = tfrgab2(x, M, L);
                    [mask,~] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect,PneiMask);
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
            RQF_tmp(Nc,it) = RQF(X(Nc,:), x_hat(Nc,:));
        end
        %RQF_tmp(:,it)
        
    end %% MC realizations
    RQF_out(:, ind_met) = mean(RQF_tmp,2);
    RQF_var_out(:, ind_met) = var(RQF_tmp,0,2);
end %% methods

   
%% Display
RQF_out
RQF_var_out