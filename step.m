
%%  Initialization for dynamics
p0 = 0; f0 = 0;                                             %initial value
% t_end = 1e-12;
% t_end1 = 1e-12;
dt = 1*1e-15;                                           %time step 
% tgrid = round(t_end/dt);                                               %number of time steps
% tgrid1 = round(t_end1/dt);
tgrid1 = 500;
tgrid = 1000;
t_end1 = tgrid1 * dt;
t_end = tgrid * dt;
p = zeros(1,round(tgrid1+1));  f = zeros(1,round(tgrid1+1));                %the function that is about to be solved
p(1) = p0; f(1) = f0;                                       %set the initial value
T = 0.1e-12;

global Omegat gamma
gamma = 0.2e12;
Omegat = zeros(1,round(tgrid1+1));
Omegat(1:round(tgrid1+1)) = sqrt(0.1)*1e12;
for n=1:tgrid1
    p(n+1) = runge_kuttap(p(n), dt, n);
    f(n+1) = runge_kuttaf(f(n), p(n), dt, n);
end
Omegat(1:tgrid+1) = 0;
for n=tgrid1:tgrid
    p(n+1) = runge_kuttap(p(n), dt, n);
    f(n+1) = runge_kuttaf(f(n), p(n), dt, n);
end

% plot((1:(501))*dt, 0.1*((1-exp(-gamma*((0:500))*dt)) / (gamma*1e-12)).^2 - abs(p(1:501)).^2);
% figure
% plot((1:(tgrid+1))*dt, abs(p).^2);
% hold on
% plot((1:(999))*dt, 0.1*((1-exp(-gamma*((0:998))*dt)) / (gamma*1e-12)).^2);
% hold on 
% plot(((501):(1001))*dt,0.1*((exp(gamma*1e-12)-1)/(gamma*1e-12)*exp(-gamma*(((500):(1000))*dt))).^2)
% figure
% plot((1:(tgrid))*dt, f);
% hold on
% plot((1:(tgrid+1))*dt, 0.2*exp(((1:(tgrid+1))*dt).^2/2))