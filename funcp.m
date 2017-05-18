function output = funcp(j, f, p)
    global gamma y Ebind A dt tstart hbar sigmat
%     output = -1i*(y'.^2.*p-1i * gamma * p/Ebind - (1-2*f).*0.1*sqrt(0.1)/Ebind/(50*dt)*hbar*(j<51))/hbar*Ebind;
  output = -1i*((y'.^2.*p- (A*2*f).*p) -3*Ebind*p - 1i * gamma * p/Ebind - (1-2*f).*((A*p)+1e-3*exp(-(j*dt+tstart)^2/(sigmat)^2)))/hbar*Ebind;
end
