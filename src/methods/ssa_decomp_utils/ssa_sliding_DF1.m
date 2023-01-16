function Y = ssa_sliding_DF1(x, N, L, nc)
% function Y = ssa_sliding_DF1(x,N,L,nc)
%
% INPUT:
% x: input signal
% N: frame length
% L: parameter of SSA
% nc: Number of components to extract
%
% Output:
% Y: the extracted components (each row vector)

N_pred = 20; %floor(N/2); %%number of sample to extrapolate


Nc = L; %nc


if (mod(N,2)==0)
 N=N+1; fprintf('warning: N modified to %d\n', N);
end

select=(N+1)/2;
Nx=length(x);
Y=zeros(Nx, Nc);

Y_pred =zeros(N_pred, Nc);

for i=1:Nx-N+1 %% for each frame  
 sx=x(i:i+N-1);

 Ytmp = zeros(N, L);
 for I = 1:L
  [Ytmp(:, I), ~] = ssa(sx, L, I); 
 end
 Nc = L;

 for ri = 1:Nc %nc  %% for each component ri = r_input, ro=r_output
  sy = Ytmp(:, ri);
  ro = ri;
  
  if (i==1)  %% begininnig
   Y(1:select,ro)=sy(1:select);
  else
   
   %% Find the best ro to connect with prediction
   scur  = sy(select+(1:N_pred)-1);  %% current values
   
   %% use prediction and match with the nearest onescur
   D = pdist2(scur', Y_pred',  'euclidean');

   %% TODO : use a threshold to set the end of a component
   [~, ro] = min(D);
   
%   plot(D)
%   pause
%    figure(33)
%    plot(scur)
%    hold on
%    plot(Y_pred(:, ro), 'r-.')
%    hold off
%    pause
   
   if (i==Nx-N+1) %% ending frame
    Y(Nx-select+1:Nx,ro) = sy(select:N);
   else %% middle frame
    Y(select+i-1, ro)=sy(select);  
   end
  end
  %% update prediction
  Y_pred(:, ro) = sy(select+(1:N_pred))';
 end
end


%% group the components

%wcorr=wCorrMat(x,L);
wcorr = pdist(Y', 'correlation' ); %'correlation'
Z = linkage(wcorr);                       
T = cluster(Z,'maxclust', nc);  

modes = zeros(Nx, nc);

for i = 1:nc        
         idx = find(T == i);
         if length(idx) > 1
          modes(:, i) = sum(Y(:, idx)');
         else
          modes(:, i) = Y(:, idx)';   
         end
end

Y = modes;



