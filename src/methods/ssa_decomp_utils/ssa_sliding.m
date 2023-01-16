function Y = ssa_sliding(x,N,L,nc)
 
if (mod(N,2)==0)
 N=N+1; fprintf('warning: N modified to %d\n', N);
end

select=(N+1)/2;
Nx=length(x);
Y=zeros(Nx,nc);
    
for i=1:Nx-N+1           
 for r=1:nc
  sx=x(i:i+N-1);
  [sy,d]=ssa(sx,L,r);            
  if (i==1)
   Y(1:select,r)=sy(1:select);
  elseif (i==Nx-N+1)
   Y(Nx-select+1:Nx,r)=sy(select:N);
  else
   Y(select+i-1,r)=sy(select);  
  end;   
 end;
end;




