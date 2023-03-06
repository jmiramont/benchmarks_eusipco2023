function [xr,if2] = nils_method(x, Ncomp, M, approach, sigma_s, clwin, return_comps,return_freq)

N = length(x);

if nargin < 3 || isempty(M)
    M  = N;
end

if nargin < 4 || isempty(approach)
    approach  = 1; % 1 = 'SR' , 2 = 'LCR'
end

if nargin < 5 || isempty(sigma_s)
    sigma_s = 0.09;
end

if nargin < 6 || isempty(clwin)
      clwin = 10;    
end

if nargin<7 || isempty(return_comps)
    return_comps = false;
end

if nargin<8 || isempty(return_freq)
    return_freq = false;
end

% -----------------------------------------------------------------------
% I added this because otherwise it uses the whole TF plane. JMM 03/02/23
if isreal(x) 
    x = hilbert(x);
end
%------------------------------------------------------------------------

[m_SR_MB, m_LCR_MB, IF_MB] = Nils_modeExtract(x, M, Ncomp,sigma_s);

if approach == 1
    X = real(m_SR_MB);
end

if approach == 2
    X = real(m_LCR_MB);
end

xr = sum(X);

if return_comps
    xr = X;
end

if2 = IF_MB/M;

if return_freq
    xr = IF_MB/M;
end