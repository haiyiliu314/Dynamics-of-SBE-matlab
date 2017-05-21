clear all
global dt t omega
n = 100;
Et = zeros(1,2*n+1); 
t = (-n:n)/(2*n+1)*10e-13*6;

dt = 1/(2*n+1)*10e-13*6;
omega = (-n:n)*2*pi*0.2413e13/(2*n+1)+2*pi*0.2413e13;
% [tgrid,omegagrid] = meshgrid(t,omega);
FFTresult = zeros(1,2*n+1);
for i=1:(2*n)
    FFTresult = runge_kuttaFT(FFTresult,i);
end
Et = exp(-((t)/(10e-13)).^2).*cos(2*pi*0.2413e13*(t));
% .*cos(2*pi*0.2413e15*(t))
plot(t,Et);
figure
% tes1 = Et.*exp(1i*omegagrid.*tgrid);
plot(omega,real(FFTresult))
