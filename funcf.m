function output = funcf(j, p)
global A dt
 output=(imag(((A*p)+0.01*exp(-(j*dt-3e-13)^2/(1e-13)^2)).* p)*2)*2*pi/6.63*10^34*4.18e-3*1.6e-19;
end

