

% test 1 :  deux sinusoides (diff amplitudes) et un chirp modulé en amplitude exponentielle

    Nx=1000;
    N=60;  L=30;   
    
    x1= real(fmconst(Nx,0.1));
    
    x2=amexpo1s(1000,1,1000).*real(fmlin(Nx,0.15,0.25));
    
    x3= real(fmconst(Nx,0.3));
    
    x=x1+2*x2+3*x3; 
    
    nc=3;  % note:  nc=4 ou 5 ne change pas le résultat.
    
   Y= ssa_sliding_JH1(x,N,L,nc);

   figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    
    % test 2   deux sinusoides de même amplitude et un chirp modulé en amplitude exponentielle
    
    x=x1+2*x2+x3; 
    
     Y= ssa_sliding_JH1(x,N,L,nc);
     
     figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    
    % test 3    deux chirps modulés (diff amplitudes) et une sinusoide
    
    x4= amexpo1s(1000,1,1000).*real(fmlin(Nx,0.3,0.35));
    
    x=x1+2*x2+3*x4;
    
    Y= ssa_sliding_JH1(x,N,L,nc); 
     
     figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
    % test 4   deux chirps modulés de même amplitude et une sinusoide
    
     x=2*x1+x2+x4;
     
      Y= ssa_sliding_JH1(x,N,L,nc); 
      
    figure
    for r=1:nc     
     subplot(nc,1,r)
     plot(Y(:,r))
     title(strcat('Y_',int2str(r)))         
    end;
    
 

 % test 5; signal audio
 
load('audio_signal.mat');
x=stmp(1:3000);
Y= ssa_sliding_JH1(x,N,L,2); 

figure
for r=1:nc     
 subplot(nc,1,r)
 plot(Y(:,r))
 title(strcat('Y_',int2str(r)))         
end;

figure  % Plot single-sided amplitude spectrum.
Fs=44100;
for r=1:nc

xx=Y(:,r);

NFFT = 2^nextpow2(length(xx)); % Next power of 2 from length of y
xf = fft(xx,NFFT)/length(xx);
fn = Fs/2*linspace(0,1,NFFT/2+1);
subplot(nc,1,r)

plot(fn,2*abs(xf(1:NFFT/2+1))) 
xlabel('Frequency (Hz)')
title(strcat('|Y_',int2str(r),'(f)|'));
end

% Note: Dans le cas du signal audio (test5), on obtient une composante sinusoidale dans Y(:,1). 
% Ce qui est intéressant dans ce cas est que la ssa classique ne permet pas
% de trouver cette composante: 
% Le spectre singulier ne donne pas d'indication sur la présence d'une
% telle composante
[y,d]=ssa(x,L,1);
figure
plot(d,'-*');

% En plus, en appliquant la ssa classique qui effectue un regroupement
% automatique par classfication hiérarchique, on ne retrouve pas cette
% composante:
Y = ssa_decomp(x, L, nc);
figure
for r=1:nc     
 subplot(nc,1,r)
 plot(Y(:,r))
 title(strcat('Y_',int2str(r)))         
end;



