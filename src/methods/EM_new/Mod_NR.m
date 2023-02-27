function wc = Mod_NR(y,wc,p,Ft,alpha)



% Perform second order gradient descent
% 
% INPUT:
% y           : Spectrogram at time n
% wc          : current mixture weight at time n
% p           : posterior distribution at time n
% Ft          : data distribution
% alpha       : Dirichlet prior parameters

% 
% OUTPUT:
% wc         : estimated mixtre weight
%
% Author: Q.Legros


%% parameters and matrix pre computation
err=10;
it=0;
[M,~,~] = size(Ft);
wc = transpose(wc);
%% Main loop
while (err>0.0001 && it<10 )
    w10=wc;       
    A = wc*Ft + (1-wc)/M;
    P = (Ft-(1/M)) ./ A;
    % Gradients
    graf = sum(sum(y.*(p'.*P)));
    H = -sum(sum(y.*(p'.*(P.^2))));

    % NR step
    wc = wc - (graf / H);
    wc = min(max(wc,1e-3),1-1e-3);
    %% Compute error
    err = sum(abs(wc-w10)./min(wc,w10));
    it = it + 1;

end
