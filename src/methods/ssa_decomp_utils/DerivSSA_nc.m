% x1=real(fmconst(60,0.1));
% x2=real(fmconst(60,0.2));
% x3=real(fmconst(60,0.3)); 
% x=x1+x2+x3; nc=3;  
% modes = DerivSSA_nc(x,L,nc); 


function modes = DerivSSA_nc(x,L,nc)

 
%[y,d]=ssa(x,L,1);

%nc=find(cumsum(d)/sum(d)>0.99,1); % estimation du nombre des composantes


N=length(x);
I1=[1 2];
dec=0;


for n=1:nc
    
    [y,d]=ssa(x,L,I1);
        I1=I1+2;
  
   if L>N/2; L=N-L; end
   K=N-L+1; 
   Tr=zeros(L,K);       
 
   for i=1:K
	Tr(1:L,i)=y(i:L+i-1);      
   end;

                  
    for i=1:N-L
    TrDr(:,i)=Tr(:,i+1)-Tr(:,i);    
    end

 % cet algorithme a un param�tre (degr� de libert�) \gamma, Here \gamma is set to 4. Otherwise it could be considered as input to the function. 

Tr =[Tr 4*TrDr]; 

S=Tr*Tr'; 
   [EigenVectors,EigenValues]=eig(S);
   %[d,i]=sort(-diag(EigenValues));  % sort(X) sorts the elements of X in ascending order.
   [d,i]=sort(diag(EigenValues),'descend');
   %d=-d;
   EigenVectors=EigenVectors(:,i); 
   
    V=(Tr')*EigenVectors;
  
   Vt=V';
   
   I=[1 2];
   
   
for r=1:nc
    
       if sum(d(I)/sum(d)>0.01)
       
   rca=EigenVectors(:,I)*Vt(I,:);
   I=I+2;
   
% Step 4: Reconstruction

   Lp=min(L,K);
   Kp=max(L,K);
   
   aux=zeros(N,1);

   for k=1:Lp-1
     for m=1:k
      aux(k)=aux(k)+rca(m,k-m+1)/k;
     
     end
   end

   for k=Lp:Kp
     for m=1:Lp
       aux(k)=aux(k)+rca(m,k-m+1)/Lp;
    
     end
   end

   for k=Kp+1:N
     for m=k+1-Kp:N+1-Kp
        aux(k)=aux(k)+rca(m,k+1-m)/(N-k+1);
 
     end
   end
   
   z(:,r+dec)=aux;
       end
end
dec=size(z,2);


end


   cor_mat=abs((corr(hilbert(z),hilbert(z))));
     
    Z = linkage(cor_mat);                       
    T = cluster(Z,'MaxClust',nc);  

  modes = zeros(N, nc);

    for i = 1:nc
        
         idx = find(T == i);
         if length(idx) > 1
         modes(:, i) = sum(z(:, idx)');
         else
          modes(:, i) = z(:, idx)';   
         end

    end
  
 [v,ind]=sort(var(modes),'descend');
 modes=modes(:,ind);
   
end

