clear all
close all
folder = './';
%% required paths 
addpath(folder);
addpath(['../signals_mat_512']);
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
SNRin = 10;

%% Apply EM method with default parameters
%best RQF for Crossing Chirps.
rng(0);
reps = 50;
M = [512 1024 2048];

L = 10;
for m =1:length(M)
    disp(M(m));
    T = round(sqrt(N));
    Pnei = floor(6*M(m)/2/pi/L);
    for i = 1:reps
        noise = randn(N,1);
        xn = sigmerge(x, noise, SNRin);
        tic()
        xr = fri_method(xn,Ncomp,M(m),L,Pnei);
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












