function [ y_hat, mask ] = EM_estim_multi(tfr, L, nb_comp, seuil, c, cl, step_Nx, ifplot  )


%% Initialize Parameters
if ~exist('c', 'var')
  c=1e-4; % TV regularization parameter
end
if ~exist('cl', 'var')
  cl=1e-2; % Laplacian spatial prior reg hyperparameter
end

if ~exist('step_Nx', 'var')
  step_Nx = 1; % Depth grid subsampling
end

% if ~exist('stepg', 'var')
%   stepg = 1e-3; % Gradient step for the M-step
% end

if ~exist('ifplot', 'var')
  ifplot = 0;  %%debug
end


%% len
N     = size(tfr,1);
M     = size(tfr,2);
Mh    = round(M/2);
  
tmp = abs(tfr(1:M/2,:)).^2;

mask  = zeros(size(tfr,1),size(tfr,2), nb_comp);
y_hat = zeros(N, nb_comp);


for i = 1:nb_comp
    
  Y = tmp;
  Ns = 1;
  [Fct]=comp_Fc(M,L);                  %% Data distribution
  Fct = repmat(Fct,[1,1,Ns]);



%% Prior choice
% reg_mu='None';
% reg_mu='Lap';
reg_mu='TV';

%% EM algorithm
tic
% [a_out,T_out]=EM_ranging(Y',Fc,Ncomp,c,si,step_Nx,stepg,nb_label);
[W_out,P,ind0]=Mod_Estim_W_EM(Y',Fct,Ns,step_Nx,reg_mu,c,cl,ifplot);
toc
  
  %% Plot reconstructions
  t_axis = (1:N)-1;
  f_axis = ((1:Mh)-1)/M;

  
%   K = 20;
%   for n = 1:size(P,2)
% 
%     [~, I] = max(P(:,n));
%     i0 = max(I-K,1);
%     i1 = min(I+K,Mh);
%     
%     if length(i0:i1) <= 0
%       continue;   
%     end
%     
%     xx = zeros(Mh,1);
%     xx(i0:i1) = P(i0:i1,n);
%     xx = [xx ; xx(end:-1:1)];
%   
%     %size(xx)
%     %size(mask(:,n,i))
%     mask(:,n,i) = xx;
%       
%   end
%   %% a modifier
  xx = zeros(size(tmp));
  xx(P.'> seuil) = 1;
  mask(:,:,i) = [xx;xx(end:-1:1,:)];

  tmp(P.'> seuil) = 0;   %% remove component

  
end

end

