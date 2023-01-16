function Y = ssa_sliding_DF0(x, N, L, nc)
% function Y = ssa_sliding_DF0(x,N,L,nc)
%
%  Use prediction to match components
%  Components are clustered before tracking
%
% INPUT:
% x: input signal
% N: frame length
% L: parameter of SSA
% nc: Number of components to extract
%
% Output:
% Y: the extracted components (each row vector)

N_pred = 10; %%number of sample to extrapolate

Nc = nc;

if (mod(N,2)==0)
 N=N+1; fprintf('warning: N modified to %d\n', N);
end

select=(N+1)/2;
Nx=length(x);
Y=zeros(Nx, Nc);

Y_pred =zeros(N_pred, Nc);

for i=1:Nx-N+1 %% for each frame  
 sx=x(i:i+N-1);
 
 Ytmp = ssa_decomp(sx, L, nc);
 

 for ri = 1:Nc %nc  %% for each component ri = r_input, ro=r_output
  sy = Ytmp(:, ri);
  ro = ri;
  
  if (i==1)  %% begininnig
   Y(1:select,ro)=sy(1:select);
  else
   
   scur  = sy(select+(1:N_pred)-1);  %% current values
   
   %% use prediction and match with the nearest onescur
   D = pdist2(scur', Y_pred');

   %% TODO : use a threshold to set the end of a component
   
   [~, ro] = min(D);
   
   if (i==Nx-N+1) %% ending frame
    Y(Nx-select+1:Nx,ro) = sy(select:N);
   
   else %% middle
    Y(select+i-1, ro)=sy(select);  
   end
   
  end
  %% update predicted samples
  Y_pred(:, ro) = sy(select+(1:N_pred))';
 end 
end



