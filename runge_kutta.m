function g = runge_kutta(g0, dt, tgrid)
g = zeros(1,tgrid+1);               %the function that is about to be solved
g(1) = g0;                        %set the initial value
for i=1:tgrid
%     t = dt * i;
    k1 = dt * func(t,g(i));
    k2 = dt * func(t + dt/2, g(i) + k1/2);
    k3 = dt * func(t + dt/2, g(i) + k2/2);
    k4 = dt * func(t + dt, g(i) + k3);
    g(i+1) = g(i) + k1/6 +k2/3 + k3/3 + k4/6;
end
plot((1:(tgrid+1))*dt, g);
hold on
plot((1:(tgrid+1))*dt, exp(((1:(tgrid+1))*dt).^2/2))

function func