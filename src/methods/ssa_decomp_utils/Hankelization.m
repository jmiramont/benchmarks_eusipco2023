function [ x_hat ] = Hankelization( X )
% [ x_hat ] = Hankelization( X )
%
% Hankelization process (cross-diagonal averaging)
%

[L,K] = size(X);

N = L + K -1;


x_hat=zeros(1,N);
Lp=min(L,K);
Kp=max(L,K);

for k=1:Lp-1
  for m=1:k
    x_hat(k)=x_hat(k)+X(m,k-m+1)/k;
  end
end

for k=Lp:Kp
  for m=1:Lp
    x_hat(k)=x_hat(k)+X(m,k-m+1)/Lp;
  end 
end 

for k=Kp+1:N
  for m=k+1-Kp:N+1-Kp
    x_hat(k)=x_hat(k)+X(m,k+1-m)/(N-k+1);
  end
end

