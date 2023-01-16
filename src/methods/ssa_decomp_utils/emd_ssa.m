clear

chemin =  './figs/emd_ssa';
if ~exist(chemin, 'dir') 
mkdir('./figs/emd_ssa/');
end


N=500;
L=40;
epsilon=0.001;

lambda0 = 0.1;
lambda1 = 0.11;     
delta_lambda=0.1;

x1=real(fmconst(N,lambda0,1));
 S(1,:)=x1;
x2= real(fmlin(N,lambda1,lambda1 + delta_lambda));

index=0;
 for logH=0.5:0.05:1   
     

 x=sigmerge(x1,x2,-20*logH);  
  x2=x-x1;
  S(2,:)=x2;

% EMD decomposition
 Y_emd = emd(x, 'fix', 50);
 Y_emd=Y_emd(1:2,:);
 [ I, s ] = match_components(S, Y_emd);

 figure(1)
  subplot(2, 2, 1);
  plot(S(1, :), 'g-');
  hold on
  plot(Y_emd(I(1),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('log_{10}(H)=%.2f\n IMF_%d, QRF=%.2f', logH, 1, s(1)));
  set(gca,'fontsize', 14); 
  xlim([0 500])
  hold off
  subplot(2, 2, 3);
  plot(S(2, :), 'g-');
  hold on
  plot(Y_emd(I(2),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('IMF_%d, QRF=%.2f', 2, s(2)));
  set(gca,'fontsize', 14); 
  xlim([0 500])
  hold off
  
% SSA decomposition 
Y_ssa = ssa_hc(x, L, 2, epsilon);
Y_ssa=Y_ssa';
[ I, s ] = match_components(S, Y_ssa);

 figure(1)
  subplot(2, 2, 2);
  plot(S(1, :), 'g-');
  hold on
  plot(Y_ssa(I(1),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('log_{10}(H)=%.2f\n SSA component_%d, QRF=%.2f', logH, 1, s(1)));
    set(gca,'fontsize', 14); 
    xlim([0 500])
  hold off
    subplot(2, 2, 4);
  plot(S(2, :), 'g-');
  hold on
  plot(Y_ssa(I(2),:), 'k-');
  set(gca,'fontsize', 14) 
  title(sprintf('SSA component_%d, QRF=%.2f', 2, s(2)));
    set(gca,'fontsize', 14); 
      xlim([0 500])
  hold off
  
index=index+1;
 saveas(gcf,sprintf('%s/emd_ssa_%d', chemin, index),'epsc');
 eps2pdf(sprintf('%s/emd_ssa_%d.eps', chemin, index));
end