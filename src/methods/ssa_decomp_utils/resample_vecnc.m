function [ vec_out ] = resample_vecnc( vec_nc, N, delta)
%[ vec_out ] = resample_vecnc( vec_nc)
%
% INPUT:
% Resample  vec_nc to c
%
%

Nframes = round(N/delta);

vec_out =  round(resample(vec_nc, N, Nframes));

end

