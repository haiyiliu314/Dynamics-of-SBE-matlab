function output = funcf(j, p)
 global Ebind dt sigmat tstart hbar

  output=(imag((1e-3*exp(-(j*dt+tstart)^2/(sigmat)^2)).*p)*2)/hbar*Ebind;
end

