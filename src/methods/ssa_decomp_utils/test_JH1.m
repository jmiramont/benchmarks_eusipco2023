

% test 1 :  deux sinusoides (diff amplitudes) et un chirp modul� en amplitude exponentielle

    Nx=1000;
    N=60;  L=30;   
    
    x1= real(fmconst(Nx,0.1));
    
    x2=amexpo1s(1000,1,1000).*real(fmlin(Nx,0.15,0.25));
    
    x3= real(fmconst(Nx,0.3));
    
    x=x1+2*x2+3*x3; 
    
    nc=3;  % note:  nc=4 ou 5 ne change pas le r�sultat.
    
   Y= ssa_sliding_JH1(x,N,L,nc);

   figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    
    % test 2   deux sinusoides de m�me amplitude et un chirp modul� en amplitude exponentielle
    
    x=x1+2*x2+x3; 
    
     Y= ssa_sliding_JH1(x,N,L,nc);
     
     figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    
    % test 3    deux chirps modul�s (diff amplitudes) et une sinusoide
    
    x4= amexpo1s(1000,1,1000).*real(fmlin(Nx,0.3,0.35));
    
    x=x1+2*x2+3*x4;
    
    Y= ssa_sliding_JH1(x,N,L,nc); 
     
     figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    % test 4   deux chirps modul�s de m�me amplitude et une sinusoide
    
     x=2*x1+x2+x4;
     
      Y= ssa_sliding_JH1(x,N,L,nc); 
      
    figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    
 % Note: je ne guarantie pas, pour N donn�, la s�paration de toutes
 %  composantes. Mais certainement en changeant le N on obtient la
 % s�paration souhait�e.
    
 


