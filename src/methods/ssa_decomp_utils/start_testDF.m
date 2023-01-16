%% Test:
%
%
%  create a multicomponent signal and recover each component using SSA
%  Compare the quality


clear all
close all


%% settings
Nx=1000;   %% length of signal
N=60;      %% windows length
L=30;

nc = 3; %% number of components


%% beginning
X = zeros(nc, Nx);
X(1,:) = real(fmconst(Nx,0.1));
X(2,:) =amexpo1s(1000,1,1000).*real(fmlin(Nx,0.15,0.25));
X(3,:) = real(fmconst(Nx,0.3));


X = X(1:nc,:);
x = sum(X);   %% signal to analyse


%% indices to consider (ignoring borders can improve reconstruction scores)
%ind = 1:Nx; %round(N/2):1:(Nx-round(N/2)+1);  


%% Method 1 classical method (not sliding)
fprintf(1, 'Classical SSA+CAH...\n');
Y1 = ssa_decomp(x, L, nc);
figure
[I1, s1] = plot_result( X, Y1, 'Classical method (not sliding SSA)');

%% Method 2 - sliding method Jinane (With tracking)
fprintf(1, 'Sliding SSA Jinane...\n');
Y2 = ssa_sliding_JH1(x,N, L, nc);
figure
[I2, s2] = plot_result( X, Y2, 'Sliding SSA Jinane 1');


%% Method 3 - classical sliding method (No tracking)
fprintf(1, 'Sliding SSA without tracking...\n');
Y3 = ssa_sliding(x,N,L,nc);
figure
[I3, s3] = plot_result( X, Y3, 'Sliding SSA (no tracking)');



%% Method 4 - Sliding method DF 0
fprintf(1, 'Sliding SSA with prediction (method 1)...\n');
Y4 = ssa_sliding_DF0_JH(x,N,L,nc);
figure
[I4, s4] = plot_result( X, Y4, 'Sliding SSA Dfourer');

%% Method 5 - Sliding method DF
fprintf(1, 'Sliding SSA with prediction (method 2)...\n');
Y5 = ssa_sliding_DF1(x,N,L,nc);
figure
[I5, s5] = plot_result( X, Y5, 'Sliding SSA Dfourer 2');











