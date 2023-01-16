function [ nc ] = nc_candidate( vec_nc )

dico = unique(vec_nc);
occ  = zeros(1, length(dico));


for i = 1:length(dico)
  occ(i) = length(find(vec_nc == dico(i)));
end

[~, idx] = max(occ);
nc = dico(idx);
end

