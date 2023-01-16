function Y = slid_ssa(x, N, delta, vec_nc, ind1, ind2, L, epsilon)
% Y = slid_ssa(x, N, delta, vec_nc, ind1, ind2, L, epsilon)
%
% INPUT:
% x: input signal
% N: frame/segment length 
% delta: the step size between successive frames
% vec_nc: vector giving the number of subcomponents to extract in each frame
% ind1, ind2: indices of the beginning (ind1) and the end (ind2) of the segment used for subcomponents matching 
                  %(the segment used for comparison between the current subcomponents and the preceding ones)  
% L: setting parameter for SSA
%
% OUTPUT:
% Y: the extracted components (each column vector)
%
%
%

if (mod(N,2)==0)
 N=N+1; fprintf('warning: N modified to %d\n', N);
end

Nx=length(x);
nb_seg = round((Nx-N)/delta)+1;
if length(vec_nc)~=nb_seg
    fprintf('warning: the vector of components number is not of suitable length')
end

if ~exist('epsilon', 'var')
 epsilon = 3*10^(-2);
end


select=(N+1)/2;
Y=zeros(Nx, max(vec_nc));
nb_comp=0;

for i=1:nb_seg% for each frame  
    
 sx = x( (i-1)*delta + (1:N) );
                    
  if nnz(sx)==0  % no signal
      continue
  end
     
Ytmp=ssa_decomp(sx, L, vec_nc(i), epsilon);

 if  i==1 ||  vec_nc(i-1) >= vec_nc(i) % we have to compare each current subcomponent to all the preceding ones
 for ri = 1:vec_nc(i)  
  
 
  sy = Ytmp(:, ri);
  ro = ri;
  
  if (i==1)  % beginning
   Y(1:select+delta-1,ro)=sy(1:select+delta-1);
  else
   
   scur  = sy(ind1 : ind2);  % current subcomponent 
   
   %% compute the euclidean distance for matching between current subcomponents and preceding ones
   D = pdist2(scur', Y_pred'); 
   
   [~, ro] = min(D);
   
   if (i==length(vec_nc)) %% ending frame
    Y(Nx-select+1:Nx,ro) = sy(select:N);
   
   else %% middle
    Y(select+  ((i-1)*delta : i*delta-1) , ro)=sy(select:select+delta-1);  
   end
   
  end
  %% update the matrix of preceding subcomponents
  Y_pred(:, ro) = sy(ind1+delta : ind2+delta)';
 end
 
 else  %  if vec_nc(i-1) < vec_nc(i) , we have to compare each preceding subcomponent to all the current ones
     
 for ri=1:vec_nc(i-1)
         
   spred = Y_pred(:, ri); % the preceding subcomponent
      
   D = pdist2(spred', Ytmp(ind1 : ind2,:)');
   
   [~, ro] = min(D);
   
   if (i==length(vec_nc)) %% ending frame
    Y(Nx-select+1:Nx,ri) = Ytmp(select:N,ro);
   
   else %% middle
    Y(select +  ((i-1)*delta : i*delta-1), ri)=Ytmp(select:select+delta-1,ro);  
   end
  Y_pred(:, ri) = Ytmp( delta+ (ind1:ind2) , ro)';
  tabro(ri)=ro; 
 end
 k=1;
 for j=1:vec_nc(i)
     if isempty(find(tabro == j)) %  the jth current subcomponent has not matched any of the preceding one 

         Y(select+  ((i-1)*delta : i*delta-1), nb_comp+k)=Ytmp(select:select+delta-1,j); % it is added to the matrix Y of components as a new component
         Y_pred(:,nb_comp+k)=Ytmp(ind1+1:ind2+1,j); % update the vector of preceding subcomponents
         k=k+1;
     end
    
 end
 
end

nb_comp=max(nb_comp,vec_nc(i));

end



