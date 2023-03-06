clear
% close all

%% required paths
folder = './';
addpath(folder);
addpath(['../signals_mat_512']);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'synchrosqueezedSTFT']));

%% Import the signals used in the paper
% load McCrossingChirps.mat
% load McSyntheticMixture5.mat
load McDampedCos.mat

N = length(x); % The signal has 1024 samples.
x = x.';
Ncomp = double(Ncomp);
X0 = comps.';

% Contaminate the signal with real Gaussian white noise.
rng(0);
SNRin = 10;

%% Apply EM method with default parameters
%best RQF for Crossing Chirps.

reps = 50;
M = [512 1024 2048];

L = 5;
for m =1:length(M)
    disp(M(m));
    T = round(sqrt(N));
    Pnei = floor(1.5*M(m)/2/pi/L); disp(Pnei);
    for i = 1:reps
        noise = randn(N,1);
        xn = sigmerge(x, noise, SNRin);        
        tic()
      % xr = brevdo_method(x, Ncomp, use_sst, Pnei, M,   L, return_comps, return_freq)
        xr = brevdo_method(xn,Ncomp,      true, Pnei, M(m),L);
        elapsed(m,i) = toc();
        qrf(m,i) = RQF(x(T+1:end-T),xr(T+1:end-T).');
    end
end

save('comp_time_variables_N_512.mat');

%%
for i = 1:3
    xx = elapsed(i,:);
    mean_elapsed(i) = mean(xx);
    mean_qrf(i) = mean(qrf(i,:)); 
end

















