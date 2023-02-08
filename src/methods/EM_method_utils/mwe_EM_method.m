clear all
close all

%% required paths
folder = './';
addpath(folder);
addpath(['signals_mat']);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

%% Import signal from file
% load McCosPlusTone.mat
% load McCrossingChirps.mat
% % load McSyntheticMixture2.mat
% % load McOnOff2.mat
% N = length(x); % The signal has 1024 samples.
% x = x.';
% Ncomp = double(Ncomp);

%% Define signal x0
N  = 500;                        %% signal length
X0(:,1) = 0.8 * fmlin(N,0.41,0.1);
X0(:,2) = fmlin(N,0.1,0.45);
x = sum(X0,2);
x = real(x); % <- Make the signal real.
Ncomp = size(X0,2);                  %% number of components

%%
% Contaminate the signal with real Gaussian white noise.
rng(0);
noise = randn(N,1);
SNRin = 30;
xn = sigmerge(x, noise, SNRin);


%% Apply EM method with default parameters
% xr = em_method(x,Ncomp,M,L,c,return_comps, return_freq)
[X] = em_method(xn,Ncomp,[],[],[],true);
xr = sum(X.',2);

%% Compute the QRF for the whole signal.
qrf = RQF(x,xr);
%qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()


%%
[H L] = roundgauss(2*N); 

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