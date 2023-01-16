clear all
close all

Ns=1000; 

a1=1; fn1=0.025;
a2=0.8; fn2=0.10; afn2=0.05;
a3=0.63; fn3=0.1625; afn3=0.0875;
Period=600;

n=(0:Ns-1)';
x1=a1*sin(2*pi*fn1*n);
x2=a2*sin(2*pi*(fn2*n+afn2*Period*sin(2*pi*n/Period)/(2*pi)));
x3=a3*sin(2*pi*(fn3*n+afn3*Period*sin(2*pi*n/Period)/(2*pi)));
sn=x1+x2+x3;

S(1,:)=x1; S(2,:)=x2; S(3,:)=x3;  


chemin =  './figs/Slidssa_Ex1/';
if ~exist(chemin, 'dir') 
mkdir('./figs/Slidssa_Ex1/');
end

 Nfft=2048; Nh=131;
 h=tftb_window(Nh,'kaiser');
 my_contrast=0.3;
 

vec_nc  = 3*ones(1, Ns);

epsilon=0.001;
delta=1;

L=40;
index=0;
for N=91:10:201;
    
ind1= (N-1)/2;
ind2= ind1+round(N/3);
Y = slid_ssa2(sn, N, delta, vec_nc, ind1, ind2, L, epsilon);

figure(9)
Y=Y';
[ I, s ] = match_components(S, Y);

for i = 1:3
  subplot(3, 1, i);
  plot(S(i, :), 'g-.');
  hold on
  plot(Y(I(i),:), 'k-');
  set(gca,'fontsize', 14) 
  if i==1
  title(sprintf('N= %d, component %d, QRF=%.2f dB', N, i, s(i)));
 % legend('reference', 'reconstruction')
  set(gca,'fontsize', 14) 
  hold off
  else
  title(sprintf('component %d, QRF=%.2f dB', i, s(i)));
    set(gca,'fontsize', 14);
 % legend('reference', 'reconstruction') 
  hold off
  end
end

%fprintf(1, 'Average RQF= %.2f dB \n', mean(s));
index=index+1;
saveas(gcf,sprintf('%s/slidssa_reconst_%d', chemin, index),'epsc');
eps2pdf(sprintf('%s/slidssa_reconst_%d.eps', chemin, index));
% 

end




