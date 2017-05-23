function output = funcp(j, f, p)
    global gamma y Ebind A dt tstart hbar sigmat
%     output = -1i*(y'.^2.*p-1i * gamma * p/Ebind- (1-2*f).*1e-3*hbar*(j<51))/hbar*Ebind;
    % - (1-2*f).*1e-3*hbar*(j<51)
  output = -1i*((y'.^2.*p- ((A*2*(f))).*p) -3*p - 1i * gamma * p/Ebind - (1-2*f).*((A*(p))+1e-3*exp(-(j*dt+tstart)^2/(sigmat)^2)))/hbar*Ebind;
end