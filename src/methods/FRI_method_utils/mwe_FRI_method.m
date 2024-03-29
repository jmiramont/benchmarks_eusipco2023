clear all
close all
folder = './';
%% required paths 
addpath(folder);
addpath(['../signals_mat']);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'FRI_lib']));

%% Import signal from file
% load McCrossingChirps.mat
% load McSyntheticMixture2.mat
load McDampedCos.mat
N = length(x); % The signal has 1024 samples0.
x = x.';
Ncomp = double(Ncomp);

% This vector tells the number of components per time sample 
% (for SSA, not used in this case).
vec_nc = double(vec_nc); 

% Contaminate the signal with real white Gaussian noise.
noise = randn(N,1);
SNRin = 5;
xn = sigmerge(x, noise, SNRin);

%% Apply the method
xr = fri_method(xn,Ncomp);

%% Compute the QRF for the whole signal.
qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()

%% Apply the method again, but recover separate components.
% xr = em_method(x,Ncomp,M,L,c,cl,step_Nx,stepg,seuil,return_comps)
xr = fri_method(xn,Ncomp,M,L,Pnei,M0,Method,return_comps, return_freq);


%%
% Show the original signal components and the recovered ones (not ordered
% by similarity).
 
[H, ~] = roundgauss(2*N); 
figure();
for i=1:Ncomp
    subplot(Ncomp+1,2,2*i-1);
    S = tfrsp(comps(i,:).',1:N,2*N,H);
    imagesc(S(1:N+1,:));
    title('Original Component: '+string(i));
    
    subplot(Ncomp+1,2,2*i);
    S = tfrsp(X(i,:).',1:N,2*N,H);
    imagesc(S(1:N+1,:));
    title('Recovered Component: '+string(i));
    
end

subplot(Ncomp+1,2,2*Ncomp+1)
S = tfrsp(sum(comps).',1:N,2*N,H);
imagesc(S(1:N+1,:));
title('Sum of Original Components: '+string(i));

subplot(Ncomp+1,2,2*Ncomp+2)
S = tfrsp(sum(X).',1:N,2*N,H);
imagesc(S(1:N+1,:));
title('Sum of Recovered Components: '+string(i));