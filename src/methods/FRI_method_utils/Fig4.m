clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.4 of the paper ESTIMATION OF OVERLAPPING MODES USING SPARSE MODELING OF SIGNALINNOVATION
%  Comparison of the RQF in the overlapping case
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
X = zeros(N,Ncomp);

X(:,1) = (fmlin(N,0.05,0.45));
X(:,2) = (fmlin(N,0.4,0.15));
X(:,3) = (fmsin(N,0.15,0.4,320,1,0.15,+1));

x0 = sum(X,2);
X = transpose(X);

[tfr] = tfrgab2(x0, M, L);       %% compute STFT
Spect = abs(tfr(1:M/2,:)).^2;



n_pad = 15; 

%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 2; % variance of the random walk in the temporal model

M0 = 10;                                    % Frequency truncation - to avoid infinite sum
F_mat = compute_F(M,L);                    % compute data distribution
F = F_mat(:,200);                          % truncate to data size

alpha = 0.5;
beta = 0.5;

%% Compute ground truth
tf0=zeros(N,1);

[tfr]  = tfrgab2(X(1,:), M, L);
tfr=(abs(tfr(1:round(M/2),:)));
for i=1:N
    [~,mm]=max(tfr(:,i));
    tf0(i,1)=mm;
end

[tfr]  = tfrgab2(X(2,:), M, L);
tfr=(abs(tfr(1:round(M/2),:)));
for i=1:N
    [~,mm]=max(tfr(:,i));
    tf0(i,2)=mm;
end

[tfr]  = tfrgab2(X(3,:), M, L);
tfr=(abs(tfr(1:round(M/2),:)));
for i=1:N
    [~,mm]=max(tfr(:,i));
    tf0(i,3)=mm;
end
tf0 = sort(tf0,2);

%% Initialization
methods_name = {'PB-\beta=0.3,Beta-div.',    'PB-\beta=0.7,Beta-div.',...
                'PB-\alpha=0.3,Renyi-div.',  'FRI',...
                'FRI_TLSD'
                };

methods_to_use = [1 2 3 4 5];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);


%% Compute MAE
RQF_out = zeros(length(SNRt), nb_methods);

for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        RQF_tmp = zeros(length(MCrep));

        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            

            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            
            switch(methods_to_use(ind_met))
                case 1  %% Beta divergence
                        beta  = 0.3; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 2;   % 2 = beta
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
                case 2  %% Beta divergence
                        beta  = 0.7; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 2;   % 2 = beta
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot);
                case 3  %% Alpha divergence
                        alpha = 0.3; % Renyi divergence hyperparameter ||  POSITIVE AND DIFFERENT TO 1
                        div   = 3;   % 2 = beta
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot);
                case 4  %% FRI
                        Method = 1;
                        [tfr] = tfrgab2(x, M, L); %% compute SST
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method);

                case 5  %% FRI TLS
                        Method = 2;
                        [tfr] = tfrgab2(x, M, L); %% compute SST
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method);
            
            end  %% switch
            [ mask ] = compMask(tf,Pnei,M/2,1);
            
            x_hat = real(rectfrgab(tfr .* mask, L, M));
            % RQF
            RQF_tmp(it) = RQF(real(x0), transpose(x_hat));

        end    %% for methods
        RQF_out(indsnr, ind_met) = mean(RQF_tmp);
    end  %% methods
 
end %% snrs
    

{'g-x' 'r-x' 'k-o' 'g-o' 'r-o'};
%% Plot -- 1
cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-o' 'g-o' 'r-o' 'k--' 'b--' 'g--' 'r--' 'k-v' 'b-v' 'g-v'};
leg = {};

figure(1)
for ind_met =  1:nb_methods
    
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, squeeze(RQF_out(:,ind_met)), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 14)
 ylabel('RQF', 'FontSize', 14)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)

grid



