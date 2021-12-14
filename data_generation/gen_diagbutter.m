%========================================================================== 
%
%   This function generate a diagonal Butterworth bandpass filter with the
%   same entries in each diagonal entry. The coefficients in numerator and
%   denominator are generated using 'butter'. The roots of B are then
%   scaled down by a small factor to keep them inside the unit circle.
%
%                   Gii(z)  = pii(z)/qii(z)
%
%                              B(1)z^n + ... + B(n-1)z + B(n)
%                           =  ------------------------------ * K
%                              A(1)z^n + ... + A(n-1)z + A(n)
%
%   After generating, a positive number K is calaculated to fix the gain of
%   the filter to 1.
%
%   INPUTS
%           N       =   filter dimension
%           n       =	order of Butterworth filter
%           band    =   lower and upper cutoff frequency in form of [a, b]
%           Fs      =   sampling frequency
%   OUTPUTS
%           Gz      =   filter's transfer function matrix
%
%	Written by ANAWAT NARTKULPAT
%
%========================================================================== 
function Gz = gen_diagbutter(N,n,band,fs)

[B,A] = butter(n,band*2/fs);
Bnew = poly(roots(B)*0.9);
Gtf = tf(Bnew,A,1/fs);
Gold = tf(B,A,1/fs);
K = abs(evalfr(Gold,exp(2j*pi*sum(band/fs))))/abs(evalfr(Gtf,exp(2j*pi*sum(band/fs))));
Gtf = K*Gtf;
%===== Plot pole-zero and bode ===
% figure; pzmap(Gold,Gtf);
% figure;  
% options = bodeoptions;
% options.FreqUnits = 'Hz';
% bode(Gold,Gtf,options);
%=================================
num = repmat({[0]},N,N); den = repmat({[1]},N,N);
Gz = tf(num,den,1/fs);

for i = 1:N
    Gz(i,i) = Gtf;
end

end

