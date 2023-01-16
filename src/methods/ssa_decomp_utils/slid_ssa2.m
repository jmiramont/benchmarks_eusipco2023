function [Y, noise] = slid_ssa2(x, N, delta, vec_nc, ind1, ind2, L, epsilon, component_dist_metric)
% Y = slid_ssa(x, N, delta, vec_nc, ind1, ind2, L)
%
% INPUTs:
% x:            input signal
% N:            frame/segment length 
% delta:        the step size between successive frames
% vec_nc:       vector giving the number of components present at each instant
% ind1, ind2:   indices of the beginning (ind1) and the end (ind2) of the segment used for subcomponents matching 
%               (should be in range [1, N])
%                 %(the segment used for comparison between the current subcomponents and the preceding ones)  
% L:            setting parameter by SSA (default is N/2)
% epsilon:      threshold used by SSA (default is 3*10^(-2))
% component_dist_metric:  choose amongs :'correlation', 'spearman', 'cityblock', 'euclidean';
%
% OUTPUT:
% Y: the extracted components (each column vector)

if ~exist('component_dist_metric', 'var')
 component_dist_metric = 'correlation'; %'correlation'; %@wcorrelation; %'spearman'; %'cityblock'; %'cityblock'; %'correlation'; %'correlation';
end

if ~exist('L', 'var')
 L = round(N/2);
end

if ~exist('epsilon', 'var')
 epsilon = 3e-2;
end

if ind2+delta > N
 warning('invalid ind2, change from %d to %d', ind2, N-delta);
 ind2 = N-delta;
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

select = (N+1)/2;                        % central sample index
Y      = zeros(Nx, max(vec_nc));         % structure to store components
noise  = zeros(1, Nx);
%Y_pred  = zeros(length(ind1:ind2), max(vec_nc));         % structure to store components

nb_comp = 0;
nbc_pred = 0;
tabro = []; %% used to store the matched components tabro(ri) = ro , where ri and ro are the components indices in 2 adjacent frames

for i = 1:nb_seg    % for each frame

  %% time indices of current frame
  I = (i-1)*delta + (1:N); I = I(I<=Nx);                   
  
  %% number of components on the previous frame
  if i > 1
   nbc_pred = nbc_cur;   
  end
  
  %% get the number of components on the current frame
  nbc_cur  = vec_nc(I(min(select, length(I))));
  %nbc_cur  = nc_candidate(vec_nc(I));
  %nbc_cur  = round(median(vec_nc(I)));
  
  sx = x(I);                                    % current frame signal
  if sum(abs(sx)) < epsilon || nbc_cur < 1    % no signal or no components
    nbc_cur = 0;
    continue;
  end
  
  %% autoSSA, TODO: to replace by a loop which reduces nbc_cur if ssa_decomp fails
  Ytmp = ssa_hc(sx, L, nbc_cur, epsilon); %% apply SSA on current frame signal
  %Ytmp = ssa_decomp(sx, L, nbc_cur, epsilon);
  %Ytmp = DerivSSA_nc(sx, L, nbc_cur);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% debug
  %size(Ytmp)
%   nbc_cur
%   plot_comp( Ytmp);
%   pause
%   plot_comp( Y);
%   pause
  
%   %% debug
%   if i > 550
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (i == 1) || (nbc_pred >= nbc_cur) %vec_nc(i-1) >= vec_nc(i) %% the number of components decreases or is the same
    for ri = 1:nbc_cur  % we have to compare each current subcomponent to all the preceding ones
      sy = Ytmp(:, ri);
      scur  = sy(ind1:ind2);  % current subcomponent used for prediction

      ro = ri;
      if i == 1  % for the first frame, directly affect each component
        Y(1:select,ro) = sy(1:select);
      else       % else, try to match current with the preceding ones     

        %% compute distance for matching between current subcomponents and preceding ones
        D = pdist2(scur.', Y_pred.', component_dist_metric);
        %D2 = pdist2(scur.', Y_pred.', 'euclidean');
        %D  = (D+D2)/2;
        %[ D ] = pdistwCorr(scur.', Y_pred.', L);
        
        [dmin, ro] = min(D);       %% ro contains the index of the closest component from scur in Y_pred

        if i == nb_seg
         I_dst = I(select):Nx;  
         I_src = select:length(sy);
         if length(I_dst) ~= length(I_src)
             m = min(length(I_dst), length(I_src));
             I_src = I_src(1:m);
             I_dst = I_dst(1:m);
         end
        else
         %I(select):(I(select)+delta-1) , ro) = sy(select:(select+delta-1));
         %I_dst = I(select):(I(select)+delta-1);
         I_src = select:(select+delta-1);
         I_dst = I(I_src);
        end
        
        if length(I_dst) ~= length(I_src)
          I_dst(1)
          I_dst(end)
          I_src(1)
          I_src(end)
          error('should never occur')   
        end
        
        %% ignore invalid matching components
%         if dmin > 0.9
%           noise(I_dst) = noise(I_dst)+sy(I_src).';
%           continue;
%         end
        
        Y(I_dst, ro) = sy(I_src);

%         if i == nb_seg  % for the last frame [OK]
% %           size(sy(select:end))
% %           size(sy(select:N))
% %           size(Y(I(select):Nx, ro))
% %           size(Y((Nx-select+1):Nx, ro))
%           Y(I(select):Nx, ro) = sy(select:end);
%           %Y((Nx-select+1):Nx, ro) = sy(select:N);
%         else                    %% middle
%           %Y(select+ ((i-1)*delta : (i*delta-1)) , ro) = sy(select:(select+delta-1));
%           Y( I(select):(I(select)+delta-1) , ro) = sy(select:(select+delta-1));
%           %Y((i-1)*delta + (1:select), ro) = sy(1:select);
%         end
      end
      %% update the matrix of preceding subcomponents
      %Y_pred(:, ro) = scur;
%       size(sy)
%       ind2+delta
%       N
      I_pred = delta + (ind1:ind2);
      I_pred = I_pred(I_pred<=length(sy));
      Y_pred(1:length(I_pred), ro) = sy(I_pred).'; %sy(delta + (ind1:ind2))';
    end
    
  else   %  if vec_nc(i-1) < vec_nc(i) ,the number of components increases
    tabro = [];
    for ri = 1:nbc_pred %size(Y_pred,2)  %nbc_pred %nbc_pred %vec_nc(i-1) %  we have to compare each preceding subcomponent to all the current ones
      spred = Y_pred(:, ri); % the preceding subcomponent 
      
      %if sum(abs(spred)) < eps
      % continue;
      %end
      
      D = pdist2(spred.', Ytmp(ind1:ind2,:).', component_dist_metric); %'correlation'
      %[ D ] = pdistwCorr(scur.', Y_pred.', L);
      %[ D ] = pdistwCorr(spred.', Ytmp(ind1:ind2,:).', L);
      
      %  D2 = pdist2(spred.', Ytmp(ind1:ind2,:).', 'euclidean');
      %  D  = (D+D2)/2;
      [~, ro] = min(D);
        
%       %% find the first element not yet affected
%       [~, Iro] = sort(D, 'descend');
%       r = 1;
%       ro = Iro(1);
%       while ~isempty(find(tabro == ro,1))
%        r = r+1;
%        ro = Iro(r);
%       end
      if i == nb_seg  %% last frame
        I_dst = Nx-select+1:Nx;
        I_src = select:N;
      else
        I_dst = I(select):(I(select)+delta-1);
        I_src = select:(select+delta-1);  
      end
      
       Y(I_dst,ri) = Ytmp(I_src,ro);

%       if i == nb_seg  %% ending frame
%         Y(Nx-select+1:Nx,ri) = Ytmp(select:N,ro);
%         %Y(((i-1)*delta + 1):Nx, ro) = Ytmp(:,ro);
%       else %% middle
%         Y(I(select):(I(select)+delta-1), ri) = Ytmp(select:(select+delta-1), ro); 
%         %Y(select +  ((i-1)*delta : i*delta-1), ri)=Ytmp(select:select+delta-1,ro); 
%         %Y((i-1)*delta + (1:select), ri) = Ytmp(1:select, ro);
%       end

      I_pred = delta + (ind1:ind2);
      I_pred = I_pred(I_pred<=length(Ytmp(:, ro)));
      
      Y_pred(1:length(I_pred), ri) = Ytmp(I_pred, ro);
      %Y_pred(:, ri) = Ytmp( delta+ (ind1:ind2) , ro)';
      %Y_pred(:, ri) = Ytmp(ind1:ind2, ro)';
      tabro(ri) = ro; 
    end
    
    %% add not yet matched components at the end
    k=1;
    for ro = 1:nbc_cur %vec_nc(i)
      if isempty(find(tabro == ro,1)) %  the jth current subcomponent has not matched any of the preceding one 
        %Y((i-1)*delta + (1:select), nb_comp+k) = Ytmp(1:select, j);
        %Y_pred(:, nb_comp+k) = Ytmp(ind1:ind2, j);
        Y(I(select):(I(select)+delta-1), nb_comp+k) = Ytmp(select:(select+delta-1), ro);  % it is added to the matrix Y of components as a new component
        Y_pred(:, nb_comp+k)                        = Ytmp(delta + (ind1:ind2), ro);      % update the vector of preceding subcomponents
        k=k+1;
      end
    end
  end
  nb_comp = max(nb_comp, nbc_cur);    %% overall number of distinct components
end %% main loop

%% remove signal where the number of components is zero
%Y(find(vec_nc == 0), :) = 0;


