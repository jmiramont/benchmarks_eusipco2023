function [ mask ] = compute_mask_comp(S, M, mu, epsilon)
% [ mask ] = compute_mask_comp(S, M, mu, epsilon)
%
%

if ~exist('mu', 'var')
    mu = 0.05;
end

if ~exist('epsilon', 'var')
    epsilon = 10^-3;
end

gw = gausswin(20);

nb_comp = size(S, 1);
N = size(S, 2);


mask = zeros(M, N, nb_comp);

for i = 1:nb_comp
 [tfr, rtfr] = tfrlmrgab(S(i,:).',1:N, M, mu);
 mask(:,:,i) = rtfr > epsilon;
 
 %% post-processing used to smooth
 tmp = conv2(mask(:,:,i), gw, 'same') / sum(gw);
 mask(:,:,i) = tmp > 0.5;
end

end
