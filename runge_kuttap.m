function g = runge_kuttap(g1, dt, i)

%     t = dt * i;
    k1 = dt * funcp(i,g1);
    k2 = dt * funcp(i, g1 + k1/2);
    k3 = dt * funcp(i, g1 + k2/2);
    k4 = dt * funcp(i+1, g1 + k3);
    g = g1 + k1/6 +k2/3 + k3/3 + k4/6;
end