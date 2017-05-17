function g = runge_kuttaFT(g1, i)
    global dt
    k1 = dt * funcFT(i);
    k2 = dt * funcFT(i+1/2);
    k3 = dt * funcFT(i+1/2);
    k4 = dt * funcFT(i+1);
    g = g1 + k1/6 +k2/3 + k3/3 + k4/6;
end

function output = funcFT(j)
    global dt freqgrid pt
    output = pt(round(j))*exp(1i*freqgrid*(dt*j));
end