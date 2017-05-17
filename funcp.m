function output = funcp(j, f, p)
    global gamma y Ebind A dt
    output = -1i*((y'.^2.*p- (A*2*f).*p) - 1i * gamma * p/Ebind - (1-2*f).*((A*p)+1e-3*exp(-(j*dt-3e-13)^2/(1e-13)^2)))*2*pi/6.63e-34*Ebind;
end
