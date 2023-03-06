clear all
close all


folder = './';
addpath('../signals_mat_512');
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));


%% Import signal from file (from the SignalBank in python).
% load McCosPlusTone.mat
% load McSyntheticMixture3.mat
load McDampedCos.mat
N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);
X0 = comps.';

% This vector tells the number of components per time sample.
vec_nc = double(vec_nc); 

% Contaminate the signal with real white Gaussian noise.
noise = randn(N,1);
SNRin = 10;
xn = sigmerge(x, noise, SNRin);

%% Apply the method: 
% xr = brevdo_method(x, Ncomp, use_sst, Pnei, M, L, return_comps, return_freq)
L = 5;
M = N;
[xr, mask, tf] = brevdo_method(xn,Ncomp,true, 15, M, L);

%% Compute the QRF
qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()


%% Apply the method again, but recover separate components.
X = brevdo_method(xn,Ncomp,true, 48, M, L, true);
[ I, s ] = match_components(comps, X);
%%
[H] = roundgauss(2*N);

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


