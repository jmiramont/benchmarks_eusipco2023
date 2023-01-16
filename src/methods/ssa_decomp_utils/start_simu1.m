clear all
close all


%%% Modele choisi   s(n) = cos(2\pi\lambda_0 n) + H \delta(n-n_0)
L = 20;


N = 500;
N_sin = 400;


s1 = zeros(1, N);

imp_pos = round((N-N_sin)/2);
imp_val = 10;
s1(imp_pos) = imp_val;


s2 = [zeros(1,N-N_sin) real(fmconst(N_sin, 0.2))'];

X = [s1;s2];


%[~,d] = ssa(s, L,1:4);

Y=ssa_decomp(s, L, 2);

subplot(121);
plot(Ytmp(:,1));
subplot(122);
plot(Ytmp(:,2));

[I4, s4] = plot_result( X, Y4, 'Sliding SSA Dfourer');
