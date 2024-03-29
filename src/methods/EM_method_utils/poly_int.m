function tf = poly_int(muc,step_r)

[N,Ns]=size(muc);
tf = zeros(size(muc));


win_close = inf*ones(size(muc));
ind = 1:size(muc,2); % indices of components
temp = repmat((1:N)',[1,Ns]);  % removed indices
for ns = 1:Ns % for loop on components
    ind2 = ind;
    ind2(ns) = [];
    for ns2 = ind2 % loop on components other than Ns
        win_close(:,ns) = abs(muc(:,ns) - muc(:,ns2)); % store abs differences
        temp(win_close(:,ns)<2*(step_r),ns) = 0;
    end
end


for ns=1:Ns
    ind_k = temp(temp(:,ns)~=0,ns); % indices of kept estimates
    ind_r = find(temp(:,ns)==0); % indices for interpolation
    if isempty(ind_k)
        tf(:,ns) = muc(:,ns);
    elseif ~isempty(ind_k)
        p = polyfit(ind_k, muc(ind_k,ns),3);
        tf(:,ns) = polyval(p,(1:N)');
        muc(ind_r,ns) = tf(ind_r,ns);
        tf(:,ns) = muc(:,ns);
    else 
        tf(:,ns) = muc(:,ns);
    end
end
tf = round(tf);