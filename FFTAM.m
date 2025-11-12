function [mod,phase,ft,faxis]=FFTAM(sig,t,T)
%% Explanation of this code
%This is a homemade FFT with normalisation and fft shift because matlab is
%apparently bad at it.

%% Validation
validateattributes(sig,      {'double'}, {'nonempty'});
validateattributes(t,        {'double'}, {'nonempty'});
validateattributes(T,       {'double'}, {'nonempty'});

%%Actual FT


L=length(t);
Fs=L/T;
%solo=transpose(hann(L)).*sig';
solo=sig';
fftest=fft(solo);
P2 = abs(fftest/L);
P2phase=angle(fftest/L);
P2full=fftest/L;
P1 = P2(1:floor(L/2)+1);
P1phase = P2phase(1:floor(L/2)+1);
P1full=P2full(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
P1phase(2:end-1) = 2*P1phase(2:end-1);
P1full(2:end-1) = 2*P1full(2:end-1);
faxis = Fs*(0:(L/2))/L;
phase=P1phase;
mod=P1;
ft=P1full;

end