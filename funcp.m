function output = funcp(j, f, p)
    global gamma A y dt Ebind
    output = -1i*(y'.^2.*p- 4*Ebind*p- (A*2*f).*p - (1-2*f).*((A*p)+1e-3*exp(-(j*dt-3e-13)^2/(1e-13)^2))- 1i * gamma * p/Ebind)*2*pi/6.63e-34*Ebind;
end
