function Y = slid_ssa0(x, N, delta, vec_nc, ind1, ind2, L, epsilon)
% Y = slid_ssa(x, N, delta, vec_nc, ind1, ind2, L)
%
% INPUTs:
% x:            input signal
% N:            frame/segment length 
% delta:        the step size between successive frames
% vec_nc:       vector giving the number of subcomponents to extract in each frame
% ind1, ind2:   indices of the beginning (ind1) and the end (ind2) of the segment used for subcomponents matching 
%               (should be in range [1, N])
%                 %(the segment used for comparison between the current subcomponents and the preceding ones)  
% L:            setting parameter by SSA (default is N/2)
% epsilon:      threshold used by SSA (default is 3*10^(-2))
%
%
% OUTPUT:
% Y: the extracted components (each column vector)

if ~exist('L', 'var')
 L = round(N/2);
end

if ~exist('epsilon', 'var')
 epsilon = 3*10^(-2);
end


if (mod(N,2)==0)
  N=N+1;
  fprintf('warning: N must be odd, modified to %d\n', N);
end

Nx=length(x);
nb_seg = round((Nx-N)/delta)+1; %% number of frame on signal

if length(vec_nc)~=Nx
  error(1, 'invalid vector of components number, should have the length of the signal')
  %warning(1,'warning: invalid vector of components number, length(vec_nc)=%d, should be %d', length(vec_nc), nb_seg);
end

select = (N+1)/2;                 % central sample index
Y=zeros(Nx, max(vec_nc));         % structure to store components

nb_comp=0;
tabro = []; %% used to store the matched components tabro(ri) = ro , where ri and ro are the components indices in 2 adjacent frames

for i=1:nb_seg   % for each frame
    
  I = (i-1)*delta + (1:N);
  I = I(I <= Nx);
  if i > 1
    nbc_pred = nbc_cur;                      % number of components on the previous frame
  end
  nbc_cur  = vec_nc((i-1)*delta + select);   % get the number of components on the current frame
  %round(median(vec_nc(I)));
  
  sx = x(I);                                    % current frame signal
  if sum(abs(sx)) < epsilon || nbc_cur == 0     % no signal or no components
    continue;
  end
 
  %% TODO: to replace by a loop which reduces nbc_cur if ssa_decomp fails
  Ytmp = ssa_decomp(sx, L, nbc_cur, epsilon); %% apply SSA on current frame signal

%   %% debug
%   if i > 129
%     i
%     
%     plot_comp( Ytmp);
%     title('Composante actuelles')
%     pause
%     
%     plot_comp( Y);
%     title('Composante enregistrees')
%     pause
%     
%     plot_comp(Y_pred);
%     title('Composante predites')
%     pause
%     nbc_cur
%     nbc_pred
%     nb_comp
%     tabro
%   end
  
  if  i == 1 ||  nbc_pred >= nbc_cur %vec_nc(i-1) >= vec_nc(i) %% the number of components decreases or is the same
    
    for ri = 1:nbc_cur  % we have to compare each current subcomponent to all the preceding ones
      sy = Ytmp(:, ri);
      ro = ri;
     
      if i == 1  % for the first frame, directly affect each component
        Y(1:select,ro)=sy(1:select);
      else       % else, try to match current with the preceding ones     
        scur  = sy(ind1 : ind2);  % current subcomponent 

        %% compute the euclidean distance for matching between current subcomponents and preceding ones
        %D = pdist2(scur', Y_pred'); 
        D = pdist2(scur', Y_pred', 'correlation');
        [~, ro] = min(D);       %% ro contains the index of the closest component from scur in Y_pred
        
        if i == nb_seg          %% for the last frame
          Y(Nx-select+1:Nx,ro) = sy(select:N);
        else %% middle
          Y(select+((i-1)*delta : i*delta-1) , ro)=sy(select:select+delta-1);  
        end
      end
      %% update the matrix of preceding subcomponents
      Y_pred(:, ro) = sy(ind1+delta : ind2+delta)';
    end
 
  else   %  if vec_nc(i-1) < vec_nc(i) ,the number of components increases
    
    for ri = 1:nbc_pred %vec_nc(i-1) %  we have to compare each preceding subcomponent to all the current ones
      spred = Y_pred(:, ri); % the preceding subcomponent 
      D = pdist2(spred', Ytmp(ind1 : ind2,:)', 'correlation');
      [~, ro] = min(D);
   
      if i == nb_seg  %% ending frame
        Y(Nx-select+1:Nx,ri) = Ytmp(select:N,ro);
      else %% middle
        Y(select +  ((i-1)*delta : i*delta-1), ri)=Ytmp(select:select+delta-1,ro);  
      end
      Y_pred(:, ri) = Ytmp( delta+ (ind1:ind2) , ro)';
      tabro(ri)=ro; 
    end
    k=1;
    for j = 1:nbc_cur %vec_nc(i)
      if isempty(find(tabro == j)) %  the jth current subcomponent has not matched any of the preceding one 
        Y(select+((i-1)*delta : i*delta-1), nb_comp+k)=Ytmp(select:select+delta-1,j); % it is added to the matrix Y of components as a new component
        Y_pred(:,nb_comp+k)=Ytmp(ind1+1:ind2+1,j); % update the vector of preceding subcomponents
        k=k+1;
      end
    end
  end
  nb_comp = max(nb_comp, nbc_cur); %vec_nc(i));
end %% main loop
