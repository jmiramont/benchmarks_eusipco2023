clear all
close all
%
%  S: matrix of input signal (line are indices of component, and columns are time indices)
%  N: length of signal
%  


%% 1 construct signal
N = 1000;
% S = zeros(2, N);
% %S(1, :) = real(fmlin(N, 0.1, 0.2));
% S(1, :) = real(fmconst(N, 0.2));
% S(2, :) = 1.5*real(fmlin(N,0.08,0.4));



a1=1; fn1=0.125;
a2=0.8; fn2=0.20; afn2=0.05;
a3=0.63; fn3=0.2625; afn3=0.0875;
Period=600;

n=(0:N)';
S(1,:)=a1*sin(2*pi*fn1*n); 
S(2,:)=a2*sin(2*pi*(fn2*n+afn2*Period*sin(2*pi*n/Period)/(2*pi)));
S(3,:)=a3*sin(2*pi*(fn3*n+afn3*Period*sin(2*pi*n/Period)/(2*pi)));

x = sum(S);

%% 2 Compute the vector containing the number of components
vec_nc = zeros(size(x));
for i = 1:size(S, 1)
  I = find(abs(S(i, :)) > eps);
  vec_nc(I) = vec_nc(I) + 1;
end

%% 3 display tfr
figure
Nh = 81; %127;% short-time window length
Nf = 256;% # of frequency bins
w = tftb_window(Nh,'Kaiser');
[sp,rs] = tfrrsp(x(:),1:length(x),Nf,w);
nfreq = (1:size(rs,1)-1)/ (2*Nf);
imagesc(1:size(rs,1),nfreq,rs(1:128,:).^0.3)
set(gca,'YDir','normal')
%set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time indices')
ylabel('frequency')
title('reassigned spectrogram of the signal')


%% 4.1 SSA settings
Nssa     = 60; %80; %71; %70;   % 81
L        = 30; %30;   %40
delta    = 1;         % step between frames (2,3,4,5 have been tested)
select   = round(Nssa/2);

% prediction and partial matching parameters
size_fen = 30;  %% 10
ind1     = select;   %select;
ind2     = ind1+size_fen-1;
epsilon = 5e-2; %5e-2; %5e-2; %3e-2; %eps; % 0.001; %3e-1; %eps

% %% 4.2 apply sliding SSA
Y1 = slid_ssa2(x, Nssa, delta, vec_nc, ind1, ind2, L, epsilon);
 %Y1 = slid_ssa2(x, Nssa, delta, vec_nc, ind1, ind2, L, epsilon);
figure
[I, s] = plot_result(S, Y1, 'Sliding SSA');

%% 4.3 apply autoSSA
Y2 = ssa_decomp(x, L, max(vec_nc), epsilon);
figure
[I, s] = plot_result(S, Y2, 'autoSSA');
