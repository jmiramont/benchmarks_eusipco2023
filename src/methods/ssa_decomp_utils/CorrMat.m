function [ cor ] = CorrMat( s, S )
% [ cor ] = CorrMat( s, S )
%
% INPUT:
% s: signal to compare of length N
% S: matrix of size nc * N containing the components to compare
%
s = s(:).';

N = length(s);
nc = size(S,1);

if N ~= size(S,2)
 error('Invalid input')   
end

cor = zeros(1, nc);

for i = 1:nc
 cor(i) = abs(sum((s-mean(s)) .* (S(i,:)-mean(S(i,:))))) / (norm(s-mean(s)) * norm(S(i,:)-mean(S(i,:))));
 %cor(i) = abs(sum(s .* S(i,:)))  / (norm(s) * norm(S(i,:)));    
end

end

