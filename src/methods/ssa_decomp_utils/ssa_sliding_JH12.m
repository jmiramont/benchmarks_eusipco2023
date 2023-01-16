

function Y = ssa_sliding_JH1(x,N,L,nc)



if mod(N,2)==0
    N=N+1;  fprintf('warning: N modified to %d\n', N);
end

    select=(N+1)/2; 
    Nx=length(x);
   
    Y=zeros(Nx,nc);  

for i=1:Nx-N+1
        
       
        sx=x(i:N+i-1);
      
        sy = DerivSSA_nc(sx,L,nc,4);  % gamma =4
        
        if abs(corr(sy(:,1),sx)) > 0.9  % si la reconstruction n'aboutit pas , je fais gamma = 0.5
            sy=  DerivSSA_nc(sx,L,nc,0.5);
        end
      if i==1
             Y(1:select,:)=sy(1:select,:);
             sy_prec=sy;        
      else
           
          
            sy_an= hilbert(real(sy));
            sy_prec_an=hilbert(real(sy_prec));
            mat_corr=abs(corr(sy_prec_an, sy_an)); % 
          
            
        for r=1:nc
                
            idx= find(mat_corr(r,:)==max(mat_corr(r,:))); 
       
         
            if idx ~= r && mat_corr(idx,idx)<0.9
               
            
             sy_prec(:,r)=sy(:,idx); 
             Y(select+i-1,r)=sy(select,idx);
             if i== Nx-N+1
               Y(Nx-select+2:Nx,r)=sy(select+1:N,idx);
             end
            
            else
                sy_prec(:,r)=sy(:,r);
                Y(select+i-1,r)=sy(select,r);
                 if i== Nx-N+1
               Y(Nx-select+2:Nx,r)=sy(select+1:N,r);
                 end
                
            end
            
        end
      end
            

   

end


