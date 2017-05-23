function g = runge_kuttaf(f1, g1, i)
    global dt
%     t = dt * i;
    k1 = dt * funcf(i,g1);
    k2 = dt * funcf(i+1/2, g1 + k1/2);
    k3 = dt * funcf(i+1/2, g1 + k2/2);
    k4 = dt * funcf(i+1, g1 + k3);
    g = f1 + k1/6 +k2/3 + k3/3 + k4/6;
end
function output = funcf(j, p)
 global Ebind dt sigmat tstart hbar
 output=(imag((1e-3*exp(-(j*dt+tstart)^2/(sigmat)^2)).*p)*2)/hbar*Ebind; 
%   output=(imag((1e-3*(j<51)).*p)*2)/hbar*Ebind;
end
