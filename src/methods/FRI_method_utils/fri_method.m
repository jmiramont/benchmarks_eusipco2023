function xr = fri_method(x,Ncomp,M,L,Pnei,M0,Method,return_comps, return_freq)
% This is the function that is called from the Python-based benchmark.

N = length(x);
%% Default Parameters
if nargin < 3 || isempty(M)
    M  = N;
end

if nargin < 4 || isempty(L)
    L  = 10; % round(N/25); % 20       %% analysis window size (in bin)
end

if nargin < 5 || isempty(Pnei)
    Pnei = 35; 
end

if nargin < 6 || isempty(M0)
    M0 = 10;  
end

if nargin < 7 || isempty(Method)
    Method = 2;
end

if nargin<8 || isempty(return_comps)
    return_comps = false;
end

if nargin<9 || isempty(return_freq)
    return_freq = false;
end


%% Apply the method
F_mat = compute_F(M,L);                    % compute data distribution
F = F_mat(:,200);                          % truncate to data size

[tfr] = tfrgab2(x, M, L); %% compute SST
Spect = abs(tfr(1:M/2,:)).^2;
[tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method);

% Get the mask and invert the tfr
[ mask ] = compMask(tf,Pnei,M/2,0);
mask(mask~=0)=1;

% Inversion of the masked STFT.
for c = 1:Ncomp
    x_hat(:,c) = real(rectfrgab(tfr .* mask(:,:,c), L, M));
    %     [xr,~] = tfristft(tfr .* mask(:,:,c),1:N,w,0);
end

if return_comps
    xr = x_hat.';
else
    [ mask ] = compMask(tf,Pnei,M/2,1);
    mask(mask~=0)=1;
    % imagesc(mask)
    xr = real(rectfrgab(tfr .* mask, L, M));
    % xr = sum(x_hat,2).'/Ncomp;
end

if return_freq
    xr = tf.'/M;
end    

    
