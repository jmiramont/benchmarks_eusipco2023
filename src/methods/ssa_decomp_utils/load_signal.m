% script load_signal
%
% generate a test signal s of length N, using variable signal = 1-6
%
% Tests signals:
%
% 1: Multicomponent signal made of 1 sinsoid and 1 linear chirps and a fmsin-chirp
% 2: Multicomponent signal made of 3 impulses
% 3: Multicomponent signal made of 3 successive sinusoids at different frequencies
% 4: Multicomponent signal made of 3 simultaneous sinusoids at different frequencies
% 5: Elementary signal made of 1 pure sinusoid at 100Hz
% 6: Real audio speech signal
%
% Recommended requirements: TFTB (http://tftb.nongnu.org/index_fr.html)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if ~exist('N', 'var')
  N = 500;    
end

if (signal == 1)
 %Nchirp  =  420; % loc_impulse1=10;loc_impulse2=50; val_impulse=15.0;  %%
 if ~exist('fmlin')
     error('Please install TFTB (http://tftb.nongnu.org/index_fr.html)')
 else
  %% Oscillating frequency  chirp
  S = zeros(3, N);
  S(1,:) = fmconst(N, 0.1);
  S(2,:) = fmlin(N,0.13,0.3);
  S(3,:) = fmsin(N,0.3,0.45,320,1,0.3,+1);
  s = sum(S);
  s = s(:);
  % 
  %  s=[zeros(N-Nchirp,1);1.0*real(s)];
  %  s(loc_impulse1)=val_impulse;
  %  s(loc_impulse2)=val_impulse;
  %  
  %% store the number of components
  %  nb_comp = zeros(size(s));
  %  nb_comp(loc_impulse1) = 1;
  %  nb_comp(loc_impulse2)       = 1;
  %  nb_comp((N-Nchirp+1):end) = 3;
   nb_comp(size(s)) = 3;
 end
elseif (signal == 2)
 N=500;
 loc_impulses= [15 150 300]; 
 
 val_impulse=10.0;
 N=1000; s=zeros(N,1); s(loc_impulses)=val_impulse;% s(150)=val_impulse; s(300)=val_impulse;
 % m=117; nfreqs(m)*Fs, figure(6); plot(t,tfr(m,:),t,rtfr(m,:)); 
 
  %% store the number of components
  nb_comp = zeros(size(s));
  nb_comp(loc_impulses) = 1;
 
elseif (signal==3)
 N=1000; L1=350; L2=300; s=[real(fmconst(L1,0.1));real(fmconst(L2,0.2));real(fmconst(N-L1-L2,0.3))];
 nb_comp = ones(size(s));

elseif (signal == 4)
 N=1000; s = real(fmconst(N, 0.1) + fmconst(N, 0.25) + fmconst(N, 0.4));
 % My_t=500; figure(6); plot(nfreqs*Fs,tfr(:,My_t),nfreqs*Fs,rtfr(:,My_t));
 nb_comp = ones(size(s)) * 3;
elseif (signal == 5)
  N=1000; s = real(fmconst(N, 0.1));
 % My_t=500; figure(6); plot(nfreqs*Fs,tfr(:,My_t),nfreqs*Fs,rtfr(:,My_t));
 nb_comp = ones(size(s));
elseif (signal == 6)
  [s, Fs] = wavread('wav/voixr.wav');
  deb = 600;
  fin = deb + Fs-1;
  s = s(deb:fin);
  
  %% override user parameters 
  N = length(s);
  
  M = 4000;
  mi = 1;mf = 700; %550;
  k = 3; L = 30;
  mu_i= 0.8; mu_f=3;    %% mu range used for LM reassignment

end