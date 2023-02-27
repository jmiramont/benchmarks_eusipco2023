function [W_out,P,ind0]=Mod_Estim_W_EM(Y,Fct,Ns,step_Nx,reg_mu,c,cl,ifplot)
    

% Main algorithm: estimate the mixture weights
% 
% INPUT:
% Y         : Spectrogram of the MCS
% Fct       : Convolution matrix of the Gaussian kernel - data distribution
% Ns        : Number of spectral component (1 since mono-componant case)
% step_Nx   : Likelihood subsampling factor (for future research, 1 for now)
% reg_mu    : Spatial regularisation (between sequential time slices)
% c         : TV spatial prior weight
% cl        : Laplacian spatial prior weight
%
% OUTPUT:
% W_out    : Mixture weight estimates
% P        : IF Posterior distribution
% ind0     : Non zeros TF instants
%
% Author: Q.Legros

% Y = Y./sum(Y,2);
%% Initialization
[N,M]=size(Y);% nb of time bins x nb of frequency bins
Nx=size(Fct,2); % nb of admissible frequencies
Nz = Nx/step_Nx; % subsampling the frequency grid for first iterations -> speed up the estimation

%% modif DF
 ind0=cell(N,1); % cell for non-empty time bins
if ifplot
    cols = {'r-.', 'g-.', 'b-.', 'c--', 'm-.', 'g-x', 'w-o'};
end


wc=0.5*ones(N,Ns); % Initialization current W
wold = wc; % Store previous iteration of W
muc=[floor(M/(Ns+1)):floor(M/(Ns+1)):M-floor(M/(Ns+1))].*ones(N,Ns); % Initialization current mu
muold = muc; % Store previous iteration of mu
sigc=floor(M/((Ns+1)*3))*ones(N,Ns); % Initialization current sigma
sigold = sigc; % Store previous iteration of sigma
    
iteEM = 0; % Number of EM iteration used to average W
pT=ones(1,Nz)./Nx; % Initialization uniform prior for mu (first EM iterations)

alpha=1.01*ones(1,Ns+1); % Initialization W priors fixed : alpha,beta (first EM iterations)


% Initialization without priors
if (strcmp(reg_mu,'TV') ||(strcmp(reg_mu,'Lap')))
    C=Mod_comp_plik(Y,Fct,wc); 
    C=C-max(C,[],2)*ones(1,Nz);
    P=exp(C);
    P=P./(sum(P,2)*ones(1,Nz)); % Normalization
    Mu=denSampling(1:Nz,P); % Gibbs sampling
    Mu = Mu.*step_Nx;
end
   
W_out=zeros(N,Ns); % Initialization output 
m_compt=1; % Initialization of EM iteration Count
Stopping_EM = 0; % EM stopping criterion
err0 = 10e9;
derr = 10e9;
errold = 10e9;
%% Main iterative process
while ((m_compt<=50 && derr>1e-3 && err0>1e-10) || m_compt<=50)
disp(['EM iteration : ',num2str(m_compt)])
    if (strcmp(reg_mu,'TV') ||(strcmp(reg_mu,'Lap'))) % TV or Lap regularization on mu
        if m_compt<3 % first iterations without regularization
           C=Mod_comp_plik(Y,Fct,wc); 
           C=C+ones(N,1)*log(pT);
           C=C-max(C,[],2)*ones(1,Nz);
           P=exp(C);
           P=P./(sum(P,2)*ones(1,Nz));
        else % Apply Lap prior
           C=Mod_comp_plik(Y,Fct,wc);  % using SEM
           [P,Mu]=compute_P_Lap_2n(C,Mu,c,cl,reg_mu,2);
        end
    else % Uniform prior
       C=Mod_comp_plik(Y,Fct,wc); 
       C=C+ones(N,1)*log(pT);
       C=C-max(C,[],2)*ones(1,Nz);
       P=exp(C);
       P=P./(sum(P,2)*ones(1,Nz));
    end

       %% Estim mu - means of the posterior
        muc=P*(1:Nz)';

    %% Maximization step
    if m_compt<3 % first iterations without priors
        wc=Mod_Mstep(Y,P,Fct,wc,alpha);
    else 
        wc=Mod_Mstep(Y,P,Fct,wc,alpha); % to update using Dirichlet prior for instance
    end

    err0 = sum(sum(abs(muc-muold)))./min(sum(sum(abs(muc))),sum(sum(abs(muold))));
    derr = norm(err0-errold);
    derr0(m_compt) = norm(err0-errold);
    wold=wc; % store previous W estimate
    errold = err0;
    
    



    %% Plots
    if ifplot
        
        figure(2)
        subplot(1,2,1)
        imagesc(Y')
        hold on
        for c = 1:Ns 
          IF = round(muc(:,c));
          h(c) = plot(1:N,IF, cols{c});
          label{c} = sprintf('mode %d', c);
        end
        legend(h, label);
        xlabel('time [s]')
        ylabel('frequency [Hz]')
        title('Data (noisy)')
        subplot(1,2,2)
        plot(1:m_compt,log(derr0))
        pause(0.01)
    
    end
    
    m_compt=m_compt+1;

    if (((m_compt>5) && (derr<1e-2)) || m_compt==20)  % Begin the last 5 iterations
        Stopping_EM = 1;
    end
    if Stopping_EM == 1
        W_out=W_out+wc;
        iteEM = iteEM + 1;
    end
end
W_out = W_out ./ iteEM; % mean of the (N_iter-N_bi) iterations


%% Plots with amplitude
SPRes = zeros(size(Y'));
for nn = 1:N
    for ns = 1:Ns
        SPRes(:,nn) = SPRes(:,nn) + W_out(nn,ns).*Fct(:,round(muc(nn,ns)),ns);
    end
end


if ifplot
%close all
figure(2)
imagesc(Y')
hold on
for c = 1:Ns 
  IF = round(muc(:,c));
  h(c) = plot(1:N,IF, cols{c});
  label{c} = sprintf('mode %d', c);
end
legend(h, label);
xlabel('time [s]')
ylabel('frequency [Hz]')
title('Final estimation')
pause(0.01)
end

