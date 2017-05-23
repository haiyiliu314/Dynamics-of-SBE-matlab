function [f, g] = runge_kuttapf(f1, g1, i)

global dt
    kg1 = dt * funcp(i,f1, g1);
    kf1 = dt * funcf(i,g1);
    kg2 = dt * funcp(i+1/2, f1 + kf1/2, g1 + kg1/2);
    kf2 = dt * funcf(i+1/2, g1 + kg1/2);
    kg3 = dt * funcp(i+1/2, f1 + kf2/2, g1 + kg2/2);
    kf3 = dt * funcf(i+1/2, g1 + kg2/2);
    kg4 = dt * funcp(i+1, f1 + kf3, g1 + kg3);
    kf4 = dt * funcf(i+1, g1 + kg3);
    g = g1 + kg1/6 +kg2/3 + kg3/3 + kg4/6;
    f = f1 + kf1/6 +kf2/3 + kf3/3 + kf4/6;
end




