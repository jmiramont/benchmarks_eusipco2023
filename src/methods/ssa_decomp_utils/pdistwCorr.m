function [ D ] = pdistwCorr(x, y, L)



x = reshape_vec(x);
y = reshape_vec(y);
[nx Nx] = size(x);
[ny Ny] = size(y);

if Nx ~= Ny
  error('x and y should have the same size')   
end

w = zeros(1,Nx);
for k = 1:Nx 
 w(k) = min([k Nx-k+1 L]);
end
w = w / sum(w);

ny = ceil(ny/2);
D = zeros(nx, ny);

for i = 1:nx
 for j = 1:ny
   D(i,j) = abs((w.* x(i,:) * y(j,:).') / (norm(w.*x(i,:)) * norm(y(i,:))));
%    if i ~= j
%     D(j,i) = D(i,j);
%    end
 end
end

D = 1-D;
end

function x = reshape_vec(x)
 if size(x,1) > size(x,2)
  x = x.';    
 end
end