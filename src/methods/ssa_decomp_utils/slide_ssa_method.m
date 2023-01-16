function xr = slide_ssa_method(x,Ncomp,Nssa,L,delta,size_win,epsilon, return_comps)
N = length(x);

if nargin<3 || isempty(Nssa)
    Nssa     = 71;
end

if nargin<4 || isempty(L)
    L        = 40;   %40
end

if nargin<5 || isempty(delta)
    delta    = 1;         % step between frames (2,3,4,5 have been tested)
end

% prediction and partial matching parameters
if nargin<6 || isempty(size_win)
    size_win = 30;  %% 10
end

if nargin<7 || isempty(epsilon)
    epsilon = 3e-1; %eps
end

if nargin<8 || isempty(return_comps)
    return_comps = false;
end


select   = round(Nssa/2);
ind1     = select;   %select;
ind2     = ind1+size_win-1;


% figure(); plot(Ncomp);

%% Apply the method
Y1 = slid_ssa2(x, Nssa, delta, Ncomp, ind1, ind2, L, epsilon);
xr = sum(Y1,2).';

if return_comps
    xr = Y1.';
end