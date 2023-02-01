clear
close all

%% required paths
folder = './';
addpath(folder);
addpath(['signals_mat']);
addpath(strcat([folder 'tools']));

%% Import signal from file
load McCosPlusTone.mat
%load McCrossingChirps.mat
% load McSyntheticMixture2.mat
%load McOnOff2.mat
N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);

% This vector tells the number of components per time sample 
% (for SSA, not used in this case).
%vec_nc = double(vec_nc); 

% Contaminate the signal with real white Gaussian noise.
noise = randn(N,1);
SNRin = 30;
xn = sigmerge(x, noise, SNRin);

%% Apply Nils method with default parameters
% xr = em_method(x,Ncomp,M,L,c,cl,step_Nx,stepg,seuil,return_comps)
%[xr] = em_method(xn,Ncomp,[],[],[],[],[],[],[],false);
M = 1024;
[  m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(xn, M, Ncomp);

%% for changing the assessed method affect X to the m_SR_Cl, m_SR_MB, m_LCR_Cl or m_LCR_MB
%% from Nils, the best one should be the spline approach (m_SR_MB or m_LCR_MB)
X = real(m_SR_MB);

xr = sum(X.',2);
% % method 1
% xr = sum(m_SR_Cl.',2);
% RQF(x,xr)
% % method 2
% xr = sum(m_SR_MB.',2);
% RQF(x,xr)
% % method 3
% xr = sum(m_LCR_Cl.',2);
% RQF(x,xr)
% % method 4
% xr = sum(m_LCR_MB.',2);
% RQF(x,xr)

%% Compute the QRF for the whole signal.
qrf = RQF(x,xr)
%qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()

%% Apply the method again, but recover separate components.
% xr = em_method(x,Ncomp,M,L,c,cl,step_Nx,stepg,seuil,return_comps)
%X = em_method(xn,Ncomp, [],[],[],[],[],[],[],true);

%%
[H L] = roundgauss(2*N); 

% Show the original signal components and the recovered ones (not ordered
% by similarity).

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