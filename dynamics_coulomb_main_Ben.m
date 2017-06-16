% clear all
global gamma A y dt freqgrid pt Ebind sigmat tstart hbar ygrid shift
%%  Initialization for dynamics
p0 = 0; f0 = 0;                                             %initial value
% t_end = 1e-12;
t_end1 = 25;  %ps
% dt1(1,j) = dt;                            %time step record
% tgrid = round(t_end/dt);                                               %number of time steps
Nt = 500000;
dt = t_end1/Nt;
ymax = 40;
N = 400;
ygrid = ymax/N;
p = zeros(N,1);  f = zeros(N,1);                %the function that is about to be solved
p(:,1) = p0; 
f(:,1) = f0;                                            %set the initial value
T = 0.1;
pt = zeros(1,Nt+1);
ft = zeros(1,Nt+1);

hbar = 4.135667662/2/pi;                                %meV*ps    
Ebind = 4.18;  %meV
omegar = Ebind/hbar;
gamma = 0.39;  %meV
time = (0:Nt)*dt;
shift = 4;

sigmat = 0.15;
tstart = -10;
Et = 1e-3*Ebind*exp(-((1:Nt+1)*dt+tstart).^2/sigmat^2);
%p(:,1) =  1e-3*Ebind*exp(-((1:N)'*dt+tstart).^2/(sigmat)^2);
Nfreq = 799; 
pfreq = zeros(1, Nfreq+1);
Efreq = zeros(1, Nfreq+1);
freqgrid =(E001'*Ebind + 4*Ebind)/hbar;          %THz
%%
%------------------------------------------------------------
N1=100;                 %number of points of phi
N2=50;                 %number of points of finer grid
A=zeros(N,N);           %Matrix for eigenvalue problem
a=zeros(1,N1);          %1:N
b=zeros(1,N1);          %discretization of phi
y=zeros(1,N);           %discretization of y
d=zeros(1,N);           %matrix of Coulumb potential
dN=ymax/N;
w=zeros(1,N1);          %weight
for i=1:N
    y(i)=dN*(i-1/2);
end
y1=y;
a=1:N1;
b=pi/(N1+1)*(a-(N1+1)/(2*pi) *sin((2*pi*a)/(N1+1)));
w=pi/(N1+1)*(1-cos(2*pi*a/(N1+1)));
d=zeros(N,N);

 for j=1:N
%%  finer grid for removing singularity
    y2=y1(j)-dN/2+((1:N2)-1/2)*dN/N2;
%%  calculate the nondiagonal terms
    for i=1:N1   %calculate the matrix as if singularity is not removed, to get nondiagonal terms
        d(j,:)=d(j,:)+(1./(sqrt(y(j)^2+y1.^2-2*y1*y(j)*cos(b(i)))) )*2.*y1*w(i);  %  integral before singularity removal
    end
    d(j,j)=0; %remove diagonal terms
%%  calculate diagonal terms
    for i1=1:N2
        d(j,j)=d(j,j)+sum((1./(sqrt(y(j)^2+y2(i1)^2-2*y2(i1)*y(j)*cos(b))))*2*y2(i1).*w)/N2;  
    end
 end
A=d;
A=A*dN/pi;
%%
%RK method to compute polarization and occupation

for n=1:Nt
    ft(n) = ygrid*y*f;
    pt(n) = ygrid*y*p;
    [f, p] = runge_kuttapf(f, p, n);
end

%%
%FFT of polarization
    ft(n+1) = ygrid*y*f;
    pt(n+1) = ygrid*y*p;

for i=1:(Nt)
    pfreq = runge_kuttaFT(pfreq,i);
end


%FFT of electric field


for i=1:(Nt+1)
    Efreq = runge_kuttaFT_E(Efreq,i);
end
FT_E = zeros(1,(Nt+1));
FT_E = sqrt(pi)*sigmat*1e-3*Ebind*exp(-1/2*sigmat^2*freqgrid.^2-tstart*1i*freqgrid);
% figure
% for i= 1:N
% plot(time,abs(p(N,:))/max(abs(p(N,:))))
% hold on
% end
% hold on
% plot(time(1:51), abs(sqrt(1e-3*Ebind/hbar*(50*dt))*((1-exp(-gamma/hbar*(0:50)*dt)) / (gamma/hbar*50*dt)))/max(abs(sqrt(1e-3*Ebind/hbar*(50*dt)*((1-exp(-gamma/hbar*(0:50)*dt)) / (gamma/hbar*50*dt)).^2))));
% hold on
% plot(time(51:end), (exp(gamma/hbar*25*dt)+exp(-gamma/hbar*25*dt))/(gamma/hbar*25*dt*50*dt)*exp(-gamma/hbar*(51:(Nt+1))*dt)/max((exp(gamma/hbar*25*dt)+exp(-gamma/hbar*25*dt))/(gamma/hbar*25*dt*50*dt)*exp(-gamma/hbar*(51:(Nt+1))*dt)))

%%
subplot(3,1,1)

plot(y, abs(p(:,end).^2))
head = sprintf('%s    |Pk|^2 at the end', date);
title(head)

subplot(3,1,2)
h = area(time, (Et/max(Et)).^2);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on

plot(time,abs(pt).^2/max(abs(pt).^2))
hold on

plot(time, abs(ft)/max(abs(ft)))
head = sprintf('%s    |P(t)|^2 ', date);
title(head)


subplot(3,1,3)
plot(time, Et)
head = sprintf('%s E(t)', date);
title(head)


figure
subplot(5,1,1)
plot(freqgrid*hbar, imag(pfreq))
title('Im[P({\omega})]')

subplot(5,1,2)
plot(freqgrid*hbar, imag(pfreq./(Efreq))/max(abs(imag(pfreq./(Efreq)))))
title( 'Im[P({\omega})/E({\omega})]')

subplot(5,1,3)
plot(freqgrid*hbar, real(Efreq))
hold on 
plot(freqgrid*hbar, real(FT_E))
title('Re[E(\omega)]')
xlabel('h\omega -(E_g+E_{1s})[meV]')

subplot(5,1,4)
plot(freqgrid*hbar, imag(Efreq))
hold on 
plot(freqgrid*hbar, imag(FT_E))
title('Im[E(\omega)]')
xlabel('h\omega -(E_g+E_{1s})[meV]')

subplot(5,1,5)
plot(freqgrid*hbar, abs(Efreq))
hold on 
plot(freqgrid*hbar, abs(FT_E))
title('|E(\omega)|')
xlabel('h\omega -(E_g+E_{1s})[meV]')