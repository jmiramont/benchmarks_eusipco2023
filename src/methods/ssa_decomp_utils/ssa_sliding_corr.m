function Y = ssa_sliding_corr(x,N,L)

% d�compose le signal x en nb_c composantes en appliquant une SSA glissante (param�tre L). La fen�tre glissante est de taille N.

% description : On applique une SSA sur un morceau du signal, que l'on fait glisser progressivement. 
% A chaque iteration i, on reconstruit � partir des premi�res nb_c paires de valeurs singuli�res, nb_c composantes.
% (Chaque composante provient  d'une paire de valeurs singuli�res)
% On compare les nouvelles composantes aux pr�c�dentes obtenues � l'it�ration i-1 en calculant les corr�lations entre celles ci. 

% <<<< la composante i est associ�e � celle qui la ressemble le plus en  terme de  corr�lation >>>>

% on retient une valeur de chaque composante, correspondant � la valeur au milieu. 

% Y est la matrice, de taille [Nx-N, nb_c], des composantes obtenues � la fin des int�rations.
% Nx est la taille du signal/vecteur x
% Y_Nx est une matrice, de taille [Nx,nb_c]. Ses colonnes contiennent les composantes obtenues mais de taille augment�e
% de Nx-N � Nx, et cela en concat�nant au d�but et � la fin de chaque  composante N/2 valeurs provenant du signal x.
 

    select=floor(N/2);
    Nx=length(x);
    sy=zeros(N,2);
    Y=zeros(Nx-N,2);


    for i=1:Nx-N+1
        
        sx=x(i:N+i-1);
        k=1;
        for r=1:2
            
            [sy(:,r),d]=ssa(sx,L,k:k+1);  
            k=k+2;
             
        end
        if i==1
             Y(i,:)=sy(select,:);
             sy_prec=sy;
        else
            
            mat_corr=abs(corr(sy_prec, sy)); % 
            idx= find(mat_corr(:,1)==max(mat_corr(:,1)));  % si pas de changements, idx=1;
  
            if length(idx)==2 || idx==2    %  conditions de changement de l'ordre des valeurs singuli�res 
                aux(:,1)=sy(:,2); aux(:,2)=sy(:,1);
                Y(i,1)=aux(select,1);
                Y(i,2)=aux(select,2);
                sy_prec=aux;                
            else
            Y(i,:)=sy(select,:); 
            sy_prec=sy;
            end
    
     
        end
    end
    
%     Y_Nx=zeros(Nx,2);
%     
%      for r=1:2
%     Y_Nx(1:select-1,r)=x(1:select-1); Y_Nx(Nx-select+1:Nx,r)=x(Nx-select+1:Nx);
%     Y_Nx(select:Nx-select,r)=Y(1:length(Y),r);
%     end
%  
   
end





