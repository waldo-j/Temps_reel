% Convolutional reverb - Offline
% T. Hueber - CNRS/GIPSA-lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% read the sample waveform
filename='./MG_list1_sent381.wav';
[x,Fs,bits] = wavread(filename);


% read the impulse response waveform
filename='./IMreverbs1/Five Columns.wav';
[imp,Fsimp,bitsimp] = wavread(filename);

% Keep only one channel 
imp_left = imp(:,1);

M = length(imp_left);
L = length(x);

step_fft = 1;

if (step_fft == 0)
  % convolution in the temporal domain (slower)
  x_conv = zeros(1,L+M-1);
  for n=1:L+M-1
    tmp = 0;
    if (n>=M)
      kmin = n-M+1;
    else
      kmin = 1;
    end
    
    if (n<L)
      kmax = n;
    else
      kmax = L;
    end
    %fprintf('kmin=%i,kmax=%i,n=%i\n',kmin,kmax,n);
    for k=kmin:kmax
      tmp = tmp + x(k)*imp_left(n-k+1);
    end
    x_conv(n)=tmp;
  end
  
  soundsc(x_conv,Fs);
else
  % FFT-base convolution (faster)
  x_conv = fconv(x,imp_left);
  
end

soundsc(x_conv,Fs);

%% END