% Example 1 - pipping bee signal analysis
clear all
close all
clc 

folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Compute_Amplitude_DF']));

Ncomp = 5;  %% number of components to extract


% load audio signal and resample to 16kHz
[s,Fs] = audioread('pipping.wav');
Fs_out = 16000;

s = sum(s,2);
s_hat = resample(s,Fs_out,Fs);

% Tfr parameter
L = 30;
M = 512;
M2 = floor(M/2);

% cropping signal
n0 = 45000;
n1 = n0+7000;
%n0+7000;  %% replace 2000 by 7000 for adressing a more complicated signal

t = ((n0:n1)-n0)/Fs_out;
f = m_axis(M)/M*Fs_out;
x = s_hat(n0:n1);

N = n1-n0+1;

% ss = s_hat(n0:n1-1);
% [tfr,stfr]  = tfrsgab2(ss, M, L);
[tfr,stfr]  = tfrsgab2(s_hat(n0:n1), M, L);
figure(1)
plot_tfr(abs(tfr(1:M2,:)).^2,t,f(1:M2))
xlabel('Time [s]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')

%% Bayesian method parameters
ds    = 2;    % variance of the random walk in the temporal model

beta  = 0.6;
alpha = 0.2;
div   = 4;                         % 1 = KL
                                   % 2 = beta
                                   % 3 = Renyi
                                   % 4 = AB-div
Pnei = 3;
PneiMask = 4;
detect = 1;



ifplot =  0;

[mask,~] = pseudoBay(tfr, Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot, detect, PneiMask);

cols = {'r-.', 'g-.', 'b-.', 'k--', 'm-.', 'g-x', 'w-o'};
figure(1)
hold on
for c = 1:Ncomp
  [ IF ] = mask2if( mask(1:M2,:,c) );
  h(c) = plot(t,IF/M*Fs_out, cols{c});
  label{c} = sprintf('mode %d', c);
end
legend(h, label);

figure()
for i = 1:Ncomp
    subplot(1,Ncomp,i)
    imagesc(mask(:,:,i));

end
