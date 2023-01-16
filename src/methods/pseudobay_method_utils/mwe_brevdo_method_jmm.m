clear all
close all


folder = './';
addpath('../signals_mat');
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));


%% Import signal from file (from the SignalBank in python).
load McCosPlusTone.mat
% load McSyntheticMixture3.mat
N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);

% This vector tells the number of components per time sample.
vec_nc = double(vec_nc); 

% Contaminate the signal with real white Gaussian noise.
noise = randn(N,1);
SNRin = 30;
xn = sigmerge(x, noise, SNRin);

%% Apply the method: 
xr = brevdo_method(xn,Ncomp,[],[],[],[],1e-8,9,false);

%% Compute the QRF
qrf = 20*log10(norm(x(100:end-100))/norm(x(100:end-100)-xr(100:end-100).'));

%% Compare recovered signal and the original (denoised) one.
figure();
plot(xr,'k','DisplayName','Recovered signal');
hold on; 
plot(x,'--g','DisplayName','Original signal'); 
legend()


%% Apply the method again, but recover separate components.
X = brevdo_method(xn,Ncomp,[],[],[],[],1e-8,9,true);

[H L] = roundgauss(2*N); 

figure();
for i=1:Ncomp
    subplot(Ncomp+1,2,2*i-1);
    S = tfrsp(comps(i,:).',1:N,2*N,H);
    imagesc(flipud(S(1:N+1,:)));
    title('Original Component: '+string(i));
    
    subplot(Ncomp+1,2,2*i);
    S = tfrsp(X(i,:).',1:N,2*N,H);
    imagesc(flipud(S(1:N+1,:)));
    title('Recovered Component: '+string(i));
end

subplot(Ncomp+1,2,2*Ncomp+1)
S = tfrsp(sum(comps).',1:N,2*N,H);
imagesc(flipud(S(1:N+1,:)));
title('Sum of Original Components: '+string(i));

subplot(Ncomp+1,2,2*Ncomp+2)
S = tfrsp(sum(X).',1:N,2*N,H);
imagesc(flipud(S(1:N+1,:)));
title('Sum of Recovered Components: '+string(i));


