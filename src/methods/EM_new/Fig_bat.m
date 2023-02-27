% Example 2 - bat signal analysis
clear all
close all

folder = './';


%% required paths 
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Modul_EM']));

Ncomp = 3;

load('bat2.mat');

N = length(x);
L = 9;

M = 512;
M2 = floor(M/2);

t = ((1:N)-1)/Fs*1000; %converted in ms
f = m_axis(M)/M*Fs;


[tfr,stfr]  = tfrsgab2(x, M, L);
spect=abs(tfr(1:M2,:)).^1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Lap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
imagesc(t,f(1:M2),spect)
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('Time [ms]','FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency [Hz]','FontSize', 12, 'FontWeight', 'bold')


ifplot = 0;
Pnei = 8;
Spect = abs(tfr(1:M/2,:)).^2;
step_r = 30;
step_v = 1;

[Fc]=comp_Fc(M,L);
Fc = Fc + eps;     %% Data distribution
[~,~,tf]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-1,step_r,step_v,ifplot,1);
[mask] = compMask(round(tf),Pnei,N,0);

cols = {'r-', 'g-', 'b-', 'k-', 'm-x', 'g-x', 'w-o'};
figure(1)
hold on
for c = 1:Ncomp
  [ IF ] = mask2if( mask(1:M2,:,c) );
  h(c) = plot(t,IF/M*Fs, cols{c});
%   h(c) = plot(t,tf(:,c), cols{c});
  
  label{c} = sprintf('mode %d', c);
    
end
legend(h, label);
axis square