clear
% close all

%% required paths
folder = './';
addpath(folder);
addpath(['../signals_mat']);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

%% Import the signals used in the paper
% load McCrossingChirps.mat
% load McSyntheticMixture5.mat
load McDampedCos.mat

N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);
X0 = comps.';

%% Define signal x0
% N  = 1024;                        %% signal length
% X0(:,1) = 0.8 * fmlin(N,0.41,0.1);
% X0(:,2) = fmlin(N,0.1,0.45);
% x = sum(X0,2);
% x = real(x); % <- Make the signal real.
% Ncomp = size(X0,2);                  %% number of components

%%
% Contaminate the signal with real Gaussian white noise.
rng(0);
noise = randn(N,1);
SNRin = 10;
xn = sigmerge(x, noise, SNRin);



%% Apply EM method with default parameters
%best RQF for Crossing Chirps.

M = N;
M2 = floor(M/2);
L = 20;

step_r= 3*ceil(sqrt(M/(pi*L)))+1; 
step_v= 4*step_r;
Pnei = floor(6*M/2/pi/L);

% Best for Crossing Chirps
% step_r= 34; step_v=78;

% Best for Synthetic Mixture
% step_r= 34; step_v=56;

tic()
[X,tf] = em_method(xn,Ncomp, M, L, [], step_r, step_v, true);
toc()
[xr,~,mask] = em_method(xn,Ncomp, M, L, Pnei, step_r, step_v);


%% Compute the QRF for the whole signal.
 [ I, s ] = match_components(comps, X);
qrf = RQF(x,xr);
%qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

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
    S = tfrsp(real(X0(:,I(i))),1:N,2*N,H);
    imagesc(S(1:N+1,:)); hold on;
    title('Original Component: '+string(i));
    
    subplot(Ncomp+1,2,2*i);
    S = tfrsp(X(i,:).',1:N,2*N,H);
    imagesc(S(1:N+1,:)); hold on;
    plot(tf(i,:)*2*N,'r');
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

figure()
S = tfrsp(sum(real(X0),2),1:N,2*N,H);
imagesc(S(1:N+1,:));
title('Sum of Original Components: '+string(i));

figure()
imagesc(mask(1:M/2+1,:));



