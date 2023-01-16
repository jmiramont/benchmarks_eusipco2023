% 
% a1=1.0; fn1=0.025;
% a2=1.3; fn2=0.10; afn2=0.05;
% a3=1.2; fn3=0.1625; afn3=0.0875;
% a4=1.1; fn4=0.275; 

clear


Nx=1000; 

a1=4; fn1=0.025;
a2=2; fn2=0.10; afn2=0.05;
a3=1; fn3=0.1625; afn3=0.0875;
a4=3; fn4=0.275; 
Period=600;

n=(0:Nx-1)';
x1=a1*sin(2*pi*fn1*n); %.*(n<600)
x2=a2*sin(2*pi*(fn2*n+afn2*Period*sin(2*pi*n/Period)/(2*pi))).*(n<700);
x3=a3*sin(2*pi*(fn3*n+afn3*Period*sin(2*pi*n/Period)/(2*pi))).*(n>200);
x4=a4*sin(2*pi*fn4*n) ; %.*(n>400);
x=x1+x2+x3+x4;

S(1,:)=x1; S(2,:)=x2; S(3,:)=x3;  S(4,:)=x4;


 Nfft=2048; Nh=131;
 h=tftb_window(Nh,'kaiser');
 my_contrast=1;
 
nb_comp = 4;
      
vec_nc  = zeros(1, Nx);
d = [1 1 201 1];
n = [Nx 699 Nx Nx];
for i = 1:nb_comp
  n_range         = d(i):n(i);
  vec_nc(n_range) = vec_nc(n_range) + 1;
end
epsilon=0.01;
delta=1;

L=40;

for N=81:10:301;
    
    
ind1= (N-1)/2;
ind2= ind1+round(N/3);
Y = slid_ssa2(x, N, delta, vec_nc, ind1, ind2, L, epsilon);


 [tfr,ti,nfreqs]=tfrsp(Y(:,1),1:Nx,Nfft,h);

figure(8);
subplot(4,1,1)
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast);
colormap(flipud(gray));
ylim([0 0.3]); grid;
set(gca,'YDir','normal','fontsize', 16) 
title(strcat('\it N=',int2str(N)),'Fontsize', 16)
 
 [tfr,ti,nfreqs]=tfrsp(Y(:,2),1:Nx,Nfft,h);

subplot(4,1,2)
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast);
colormap(flipud(gray));
ylim([0 0.3]); grid;
set(gca,'YDir','normal','fontsize', 16) 

 [tfr,ti,nfreqs]=tfrsp(Y(:,3),1:Nx,Nfft,h);

subplot(4,1,3)
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast);
colormap(flipud(gray));
ylim([0 0.3]); grid;
set(gca,'YDir','normal','fontsize', 16) 

  [tfr,ti,nfreqs]=tfrsp(Y(:,4),1:Nx,Nfft,h);

subplot(4,1,4)
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast);
colormap(flipud(gray));
ylim([0 0.3]); grid;
set(gca,'YDir','normal','fontsize', 16) 
pause(3)

figure(9)
[I, s] = plot_result(S, Y.');
fprintf(1, 'Average RQF= %.2f dB \n', mean(s));
pause(3)

end


