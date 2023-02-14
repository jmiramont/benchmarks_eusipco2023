clear all
close all

%%% script for searching the best set of parameters
%step_r / step-v

nval1 = 10;
nval2 = 10;

v1_min = 1;
v1_max = 100;

v2_min = 1;
v2_max = 100;


rqf_res = ones(nval1,nval2) * nan;


%% required paths
folder = './';
addpath(folder);
addpath(['signals_mat']);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

%% Import signal from file
load McCrossingChirps.mat
% load McSyntheticMixture5.mat
% load McDampedCos.mat

N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);
X0 = comps.';

%% Define signal x0
% N  = 500;                        %% signal length
% X0(:,1) = 0.8 * fmlin(N,0.41,0.1);
% X0(:,2) = fmlin(N,0.1,0.45);
% x = sum(X0,2);
% x = real(x); % <- Make the signal real.
% Ncomp = size(X0,2);                  %% number of components

%%
% Contaminate the signal with real Gaussian white noise.
rng(0);
noise = randn(N,1);
SNRin = 20;
xn = sigmerge(x, noise, SNRin);

v1_rng = linspace(v1_min,v1_max, nval1);
v2_rng = linspace(v2_min,v2_max, nval2);

best_v = [-1 -1];
best_qrf = -inf;
qrf = zeros([nval1,nval2]);
for i = 1:nval1
    parfor j = 1:nval2
      step_r = round(v1_rng(i));
      step_v = round(v2_rng(j));
      %% Apply EM method with default parameters
      %20 / 60 works
      [X] = em_method(xn,Ncomp,[],[],[],step_r, step_v, true);
      xr  = sum(X.',2);

      %% Compute the QRF for the whole signal.
      qrf(i,j) = RQF(x,xr);
      qrf(i,j)
      
%       if qrf > best_qrf
%         best_v = [step_r step_v]
%         best_qrf = qrf
%         
%         rqf_res(i,j) = qrf;
%       end
    end
end

figure
imagesc(v1_rng,v2_rng, rqf_res)
hold on
plot(best_v(1),best_v(2), 'kx')
xlabel('step_r')
ylabel('step_v')

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()


%%
[H, L] = roundgauss(2*N); 

% Show the original signal components and the recovered ones (not ordered
% by similarity).

figure();
for i=1:Ncomp
    subplot(Ncomp+1,2,2*i-1);
    S = tfrsp(real(X0(:,i)),1:N,2*N,H);
    imagesc(S(1:N+1,:));
    title('Original Component: '+string(i));
    
    subplot(Ncomp+1,2,2*i);
    S = tfrsp(X(i,:).',1:N,2*N,H);
    imagesc(S(1:N+1,:));
    title('Recovered Component: '+string(i));
end

subplot(Ncomp+1,2,2*Ncomp+1)
S = tfrsp(sum(real(X0),2),1:N,2*N,H);
imagesc(S(1:N+1,:));
title('Sum of Original Components: '+string(i));

subplot(Ncomp+1,2,2*Ncomp+2)
S = tfrsp(sum(X).',1:N,2*N,H);
imagesc(S(1:N+1,:));
title('Sum of Recovered Components: '+string(i));