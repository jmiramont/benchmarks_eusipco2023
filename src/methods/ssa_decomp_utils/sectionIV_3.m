clear

Ns=1000;
spread=0.45;

chemin =  './figs/';
if ~exist(chemin, 'dir') 
mkdir('./figs/');
end


%% 
nbc = 2;
S       = zeros(nbc, Ns);        
vec_nc  = zeros(1, Ns);
d = [20 220] ;
q = [620 720];
for i = 1:nbc
  n_range         = d(i)+(1:q(i));
  vec_nc(n_range) = vec_nc(n_range) + 1;
  S(i, n_range)   = fmlin(q(i),0.1,0.4).*amgauss(q(i),q(i)/2,spread*q(i));
end
S = real(S);

sn = sum(S).';


%% 2: signal plot

figure(1);
subplot(2,1,1)
plot(S(1,:))
set(gca,'fontsize', 16)
xlabel('time index')
ylabel('\it x^{(1)}')
subplot(2,1,2)
plot(S(2,:))
 set(gca,'fontsize', 16)
xlabel('time index')
ylabel('\it x^{(2)}')
fname=strcat(chemin, 'signal2'); 
saveas(gcf,fname,'epsc'); 
eps2pdf(strcat(fname,'.eps'));

 
%[tfr,ti,nfreqs]=tfrsp(sn,1:Ns,Nfft,h);
% figure
% imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast);
% colormap(flipud(gray));
% %ylim([0 0.5]);
% grid;
% set(gca,'YDir','normal','fontsize', 16) 
% title('Spectrogram of \it s_n')
% ylabel('normalized frequency')
% xlabel('time index')


%% 3: sliding SSA - mode extraction
chemin1 =  './figs/Slidssa_Ex2/';
if ~exist(chemin1, 'dir') 
mkdir('./figs/Slidssa_Ex2/');
end

L=40;
epsilon = 0.001;
index=0;
for N=81:10:151
    
select   = round(N/2);
ind1= (N-1)/2;
ind2= ind1+round(N/3);
delta=1;

Y = slid_ssa2(sn, N, delta, vec_nc, ind1, ind2, L, epsilon);

 figure(2)
 Y=Y';
[ I, s ] = match_components(S, Y);
subplot(2, 1, 1);
  plot(S(1, :), 'g-.');
  hold on
  plot(Y(I(1),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('N= %d, component %d, QRF=%.2f', N, 1, s(1)));
 hold off
subplot(2, 1, 2);
  plot(S(2, :), 'g-.');     
  hold on
  plot(Y(I(2),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('component %d, QRF=%.2f', 2, s(2)));
   hold off 
fprintf(1, 'Average RQF= %.2f dB \n', mean(s));
index=index+1;
saveas(gcf,sprintf('%s/slidssa_reconst_%d', chemin1, index),'epsc');
eps2pdf(sprintf('%s/slidssa_reconst_%d.eps', chemin1, index));

end

% SSA decomposition 
chemin2 =  './figs/videoEx2/';
if ~exist(chemin2, 'dir') 
mkdir('./figs/videoEx2/');
end
index=0;
 
Nfft = 100;
h = tftb_window(50,'Kaiser');
for L=40:10:200;   
    
Y= ssa_hc(sn, L, nbc, epsilon);
 
figure(3)
subplot(2,2,1)
plot(Y(:,1))
set(gca,'fontsize', 14) 
ylabel('\it y^{(1)}')
 title(sprintf('L = %d', L));

subplot(2,2,2)
[tfr,ti,nfreqs]=tfrsp(Y(:,1), 1:Ns , Nfft ,h);
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast)
colormap(flipud(gray));
ylim([0 0.5]); grid;
set(gca,'YDir','normal','fontsize', 14) 
xlim([0 1000]);
title('Spectrogram')


subplot(2,2,3)
plot(Y(:,2))
set(gca,'fontsize', 12) 
ylabel('\it y^{(2)}')

subplot(2,2,4)
[tfr,ti,nfreqs]=tfrsp(Y(:,2),1:Ns , Nfft ,h);
imagesc(ti, nfreqs(1:Nfft/2) , tfr(1:Nfft/2,:).^my_contrast)
colormap(flipud(gray));
ylim([0 0.5]); grid;
set(gca,'YDir','normal','fontsize', 14) 
xlim([0 1000]);


index=index+1;
 saveas(gcf,sprintf('%s/ssa_reconst_%d', chemin2, index),'epsc');
 eps2pdf(sprintf('%s/ssa_reconst_%d.eps', chemin2, index));


end
