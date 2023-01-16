clear all
close all


%% signal length
signal = 3;
load_signal_ssa;
T = 1:length(s);

%% Compute the vector containing the number of components
vec_nc = zeros(size(s));
for i = 1:size(S, 1)
  I = find(abs(S(i, :)) > eps);
  vec_nc(I) = vec_nc(I) + 1;
end




%% 1 SSA
Nssa     = 70;  % 81
L        = 30;   %40
delta    = 1;         % step between frames (2,3,4,5 have been tested)
select   = Nssa/2;

%% range used for prediction and partial matching
size_fen = 20;  %% 10
ind1     = select;
ind2     = ind1+size_fen-1;

epsilon = 3e-2;

%% 1 Sliding SSA
Y = slid_ssa2(s, Nssa, delta, vec_nc, ind1, ind2, L, epsilon);
[I1, s1] = plot_result( S, Y, 'Sliding-SSA');

%% 1 Classical SSA
Y_ssa = ssa_hc(s, L, max(vec_nc), epsilon);  %  (s, Nssa, delta, vec_nc, ind1, ind2, L, epsilon);
figure
[I2, s2] = plot_result( S, Y_ssa, 'Classical-SSA');


%% 2 EMD
Y_emd = emd(s, 'MAXMODES', max(1,max(vec_nc)-1));
figure
[I3, s3] = plot_result( S, Y_emd, 'EMD');


%% 3 Synchrosqueezing





% 
% 
% %% 1) time-frequency distribution
% Nf = 256;           % # of frequency bins
% Nh = 81;            %127;% short-time window length
% w = tftb_window(Nh,'Kaiser');
% %w = tftb_window(Nh,'Gauss');
% 
% [sp,rs] = tfrrsp(s,T,Nf,w);
% figure(1)
% imagesc(flipud(rs(1:128,:)))
% set(gca,'YTick',[]);set(gca,'XTick',[])
% xlabel('time')
% ylabel('frequency')
% title('signal')
