% LigoSignal analysis using SSA
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)


clear all
close all


selected_signal = 1;  %%choose the signal to process: can be 1 or 2

if selected_signal == 1
    load fig1-observed-Livingston.txt
    t1=fig1_observed_Livingston(:,1);
    sig1=fig1_observed_Livingston(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal1';
elseif selected_signal == 2

    load fig1-observed-Hanford.txt
    t1=fig1_observed_Hanford(:,1);
    sig1=fig1_observed_Hanford(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal2';
else
    error('Unknown signal')
end


s = sig1;

%% 1 SSA settings
Nssa     = 300;  % 81
L        = 100;   %40
delta    = 1;         % step between frames (2,3,4,5 have been tested)
select   = Nssa/2;

%% range used for prediction and partial matching
size_fen = 20;  %% 10
ind1     = select;
ind2     = ind1+size_fen-1;

epsilon = 0.15; %3e-1; %% should be taken in [0, 1]


%% tfr representation
figure
Nh = 301; %81; %127;% short-time window length
Nf = 2048; %256;% # of frequency bins
w = tftb_window(Nh,'Kaiser');
[sp,rs] = tfrrsp(s,1:length(s),Nf,w);
imagesc(flipud(rs(1:60,:).^0.2)) %1:128  % 60/2048 * Fe1 = 480Hz
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('reassigned spectrogram of the signal')

nb_comp = 2; %% 2 components to extract

vec_nc = ones(1,length(s))*nb_comp;

%% ssa classique
Y = ssa_decomp(s, L, nb_comp, epsilon);
figure
plot_comp(Y)
subplot(nb_comp+1,1,1)
title('Classical SSA')

figure
plot(t1, sig1, 'g-.')
hold on
plot(t1, Y(:,2), 'k')
xlabel('time (s)', 'FontSize', 16)
legend('reference signal', 'reconstructed signal');
title('Resulting signal using classical SSA')


%% sliding-ssa
Y2 = slid_ssa(s, Nssa, delta, vec_nc, ind1, ind2, L, epsilon);
figure
plot_comp(Y2)
subplot(nb_comp+1,1,1)
title('Sliding-SSA')

figure
plot(t1, sig1, 'g-.')
hold on
plot(t1, Y2(:,2), 'k')
xlabel('time (s)', 'FontSize', 16)
legend('reference signal', 'reconstructed signal');
title('Resulting signal using sliding-SSA')



