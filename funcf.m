function output = funcf(t, p)
global Omegat
omega_g = 3.2 * 1.6e-19 / (6.63e-34);
% gamma = 0.2e12;
% output=2*imag(Omegat(1) * exp(1i * omega_g * t) * p)/(6.63e-34);
output=imag(Omegat(t) * p)*2;
end