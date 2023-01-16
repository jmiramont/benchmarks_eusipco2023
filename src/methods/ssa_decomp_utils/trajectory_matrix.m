function [ X, K, N ] = trajectory_matrix( x, L )
% [ X, K, N ] = trajectory_matrix( x, L )
%
%
%

N=length(x);

if L>N/2
    L=N-L;
end

K=N-L+1;

X=zeros(L,K);                 % [L,K]=size(X)

for i=1:K
    X(1:L,i)=x(i:L+i-1); 
end;
   
end

