clear all
close all


L=22;
N = 100;
H = 7;
phi = pi/6;

lambda = 0.3;
w = 2*pi*lambda;

n = 1:N;
n0 = L + 22;

%% input signals x:sinusoid , y: impulse
x = exp(1j * w * n + 1j*phi); %real(exp(1j * w * n + 1j*phi));  %
y = zeros(1, N);y(n0) = H;
z = x+y;
%x = real(x);
%s = real(x + y);
%s = x + y;
%S = ssa_decomp(s, L, 2);

%% trajectory matrices of each part
[X, K] = trajectory_matrix(x, L);
[Y] = trajectory_matrix(y, L);
[Z] = trajectory_matrix(z, L);

%% eigenvectors (U and V from SSA)
U1 = exp(1j * w * (0:(L-1))).';
V1 = exp(1j * w * (0:(K-1)))';

%% the eigenvalue (Rank(X) = 1)
S1 = exp(1j*w+1j*phi);

%% reconstructed trajectory matrix of S1
X_hat = S1 * U1 * V1';

%% 
Emat = U1 * (Y' * U1)';

E = Hankelization( Emat )/H

x_hat = Hankelization( X_hat );

alpha = zeros(1, N);
alpha(n0) = 1;
for i = 1:(n0-1)
 alpha(n0-i) = alpha(n0-i+1)  - 1/L;
 alpha(n0+i) = alpha(n0+i-1) - 1/L;
end
alpha(alpha<= 0) = 0;

plot(abs(E))
hold on
plot(real(E), 'k-')
plot(alpha, 'r-.')
xlabel('sample index')
legend('real(x_n)', '|x_n|', '\alpha_i')


z     = cos(w * n + phi) + real(E); %2 * real(E)/norm(E)^2;


s     = real(x+y);
[z_hat d] = ssa(s,L,1:2);

figure
plot(z)
hold on
plot(z_hat, 'r-.')
