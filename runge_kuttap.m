function g = runge_kuttap(f1, g1, i)
    global dt
%     t = dt * i;
    k1 = dt * funcp(i,f1, g1);
    k2 = dt * funcp(i+1/2, f1, g1 + k1/2);
    k3 = dt * funcp(i+1/2, f1, g1 + k2/2);
    k4 = dt * funcp(i+1, f1, g1 + k3);
    g = g1 + k1/6 +k2/3 + k3/3 + k4/6;
end
function output = funcp(j, f, p)
    global gamma y Ebind A dt tstart hbar sigmat
%     output = -1i*(y'.^2.*p-1i * gamma * p/Ebind - (1-2*f).*0.1*sqrt(0.1)/Ebind/(50*dt)*hbar*(j<51))/hbar*Ebind;
  output = -1i*((y'.^2.*p- (A*2*f).*p) -3*Ebind*p - 1i * gamma * p/Ebind - (1-2*f).*((A*p)+1e-3*exp(-(j*dt+tstart)^2/(sigmat)^2)))/hbar*Ebind;
end