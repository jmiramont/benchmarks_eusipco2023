function C=Mod_comp_plik(Y,Fc,wc)


% Compute the log likelihood 
% 
% INPUT:
% Y        : Observation
% Fc       : Convolution matrix of the Gaussian kernel - data distribution
% wc       : current mixture weights

% OUTPUT:
% C        : Log likelihood
%
% Author: Q.Legros


[N,M]=size(Y);
[~,Nx,~]=size(Fc);
C=zeros(N,Nx);
% for n=1:N
%     if sum(Y(n,:))~=0   
%         A =(1-sum(wc(n,:)))/M;
%         for ns = 1:Ns
%             A = A+(wc(n,ns)*Fc(:,:,ns));
%         end
%         C(n,:)=sum((Y(n,:)'*ones(1,Nx)).*log(A),1);
%     end
% end
for n=1:N
    A = wc(n)*Fc+(1-wc(n))/M;
    tmp = Y(n,:)'*ones(1,Nx);
    C(n,:)=sum(tmp.*log(A),1);
end
