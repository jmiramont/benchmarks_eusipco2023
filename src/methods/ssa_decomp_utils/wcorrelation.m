function [ D ] = wcorrelation( X, Y )

X = X(:).';
Y = Y(:).';
N = length(X);
% w = zeros(1,N);
% for k = 1:Nx 
%  w(k) = min([k N-k+1 L]);
% end
w = hanning(N).';
w = w / sum(w);

denum = (norm(w.*X) * norm(Y));
if abs(denum) < eps
 D = 0;   
else
 D = abs((w.* X * Y.') / denum);
end

D = 1-D;
end

