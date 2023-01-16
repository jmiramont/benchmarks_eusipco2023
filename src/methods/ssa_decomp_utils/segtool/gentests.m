function [s s1 s2 s3] = gentests()
%% gentests : generates a multicomponent signal.
% Input : No
% Outputs : s=s1+s2+s3, a 3-components signal


% Parameters
N = 1024;
T = 1;
t = linspace(0,T-T/N,N);
un = ones(size(t));



%% Multicomponent signal
a1 = 2*un+1.3*cos(4*pi*t);
phi1 = 2*pi*77;
f1=30*sin(3*pi*t);
fi1 = phi1+30*3*pi*cos(3*pi*t);

a2 = 1*un+0.5*sin(3*pi*t+0.2*un);
phi2 = 2*pi*46;
f2=21*sin(3*pi*t);
fi2 = phi2+21*3*pi*cos(3*pi*t);

a3 = 2*exp(-22*(t-T/2*un).^2);
phi3 = 2*pi*8*un;
a3 = 2*exp(-22*(t-T/2*un).^2);
phi3 = 2*pi*6*un;

a1=1;a2=1;a3=1;
s1 = a1.*sin(3*(phi1*t+f1));
s2 = a2.*sin(3*(phi2*t+f2));
s3 = a3.*sin(3*(phi3.*t));
s = s1+s2+s3;



% % Figures : plots each component
% indplot = 1:12:N;
% figure();
% hold on;
% scatter(t(indplot),fi1(indplot)/(2*pi),10,a1(indplot));
% scatter(t(indplot),fi2(indplot)/(2*pi),10,a2(indplot));
% scatter(t(indplot),phi3(indplot)/(2*pi),10,a3(indplot));
% xlabel('time');
% ylabel('frequency');
% colormap(flipud((gray(128))));
% hc = colorbar('XTick',1,'XTickLabel','Amplitude');
% title('Test signal s');
% 
% % Figures : TF representation
% figure();
% subplot(4,1,1);plot(t,s);xlabel('t');ylabel('s');
% subplot(4,1,2);plot(t,s1);xlabel('t');ylabel('s_1');
% subplot(4,1,3);plot(t,s2);xlabel('t');ylabel('s_2');
% subplot(4,1,4);plot(t,s3);xlabel('t');ylabel('s_3');
% % 




end
