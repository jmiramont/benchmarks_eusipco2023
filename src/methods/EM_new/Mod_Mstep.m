function wc=Mod_Mstep(Y,P,Fc,wc,alpha)


% Perform M-step
% 
% INPUT:
% Y           : Observation
% P           : Posterior distribution
% Fc          : data distribution
% wc          : current mixture weights
% alpha       : Dirichlet prior hyperparameter

% 
% OUTPUT:
% wc         : estimated mixture weights
%
% Author: Q.Legros

[M,Nx,Ns]=size(Fc);
[N,~]=size(Y);
for n=1:N % Update pixelwise
    if sum(Y(n,:))==0 % If spectrogram is zeros value : max the posterior relies to max the prior
        wc(n,:)=(alpha(1:Ns)-1)/(sum(alpha)-length(alpha)); % mode of the Dirichlet distribution
    else % else : max the posterior relies to max the prior
        
        
    ind=find(P(n,:)>0.1*max(P(n,:))); % Select useful value of the posterior - speed-up the estimation
    p=P(n,ind)';

    wc(n,:) = Mod_NR(Y(n,:)',wc(n,:),p,Fc(:,ind,:),alpha);

    end
end