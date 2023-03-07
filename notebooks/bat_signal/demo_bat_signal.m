% Example of methods with bat signal 
% Bat chirp signal: Digitized 2.5 microsecond echolocation pulse emitted 
% by the Large Brown Bat, Eptesicus Fuscus. There are 400 samples;
% the sampling period was 7 microseconds. 
% -------------------------------------------------------------------------
clear all; close all;

% Load the signal
load batsig.mat
%fs = 1/7e-6;
% x = batsig-mean(batsig);
x = batsig;
Ncomp = 3;
%% Spectrogram:
N = length(x);
Nfft = 2*N;
H = roundgauss(2*N); 
S = tfrstft(x,1:N,2*N,H);
figure();
imagesc(flipud(abs(S(1:N+1,:))));
% colormap jet

%% EM method
folder = './';
addpath(folder);
addpath('../../src/methods/EM_method_utils')
step_r = 30;
step_v = 1;
[~,~,mask] = em_method(x,Ncomp,[],[],15,step_r,step_v,[],[],0);
mask_EM = mask(1:N/2,:);
figure()
imagesc(mask_EM);
save('mask_EM.mat','mask_EM')

%% PB method
addpath('../../src/methods/PB_method_utils/')
[tf,mask] = pb_method(x, Ncomp, false, [], 0.7, 0.3, [], 15, [], [], [], [],true);
mask_PB = mask(1:N/2,:);
figure()
imagesc(flipud(mask_PB)); hold on;
plot((0.5-tf(1,:))*N,'-r');
plot((0.5-tf(2,:))*N,'-r');
plot((0.5-tf(3,:))*N,'-r');
save('mask_PB.mat','mask_PB')

