function output = funcp(t, p)
global Omegat gamma
omega_g = 3.2 * 1.6e-19 ;
% output=(omega_g - 1i * gamma) * p -Omegat(1);
output=-1i*(( - 1i * gamma) * p -Omegat(t));
end
