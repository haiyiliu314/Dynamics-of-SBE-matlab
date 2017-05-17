clear all
global gamma A y dt freqgrid pt Ebind
%%  Initialization for dynamics
p0 = 0; f0 = 0;                                             %initial value
% t_end = 1e-12;
t_end1 = 30e-13;  
% dt1(1,j) = dt;                            %time step record
% tgrid = round(t_end/dt);                                               %number of time steps
Nt = 200;
dt = t_end1/Nt;
ymax = 4;
N = 1000;
ygrid = ymax/N;
p = zeros(N,round(Nt+1));  f = zeros(N,round(Nt+1));                %the function that is about to be solved
p(:,1) = p0; f(:,1) = f0;                                       %set the initial value
T = 0.1e-12;
pt = zeros(1,Nt+1);
hbar = 6.624e-34/2/pi;
Ebind = 4.18e-3*1.6e-19;
omegar = Ebind/hbar;
gamma = 0.39e-3*1.6e-19;
time = (0:Nt)*dt;
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
A=-1*A*dN/pi;
%%
%RK method to compute polarization and occupation

for n=1:Nt
    p(:,n+1) = runge_kuttap(f(:,n), p(:,n), n);
    f(:,n+1) = runge_kuttaf(f(:,n), p(:,n), n);
end

%%
%FFT of polarization
for j = 1:(Nt+1)
    pt(j) = ygrid*y*p(:,j);
end
pfreq = zeros(1,Nt+1);
freqgrid = ((0:Nt)/(Nt+1)* 2*omegar-omegar)*16;
for i=1:(Nt)
    pfreq = runge_kuttaFT(pfreq,i);
end


%FFT of electric field
Et = exp(-((1:Nt+1)*dt-3e-13).^2/(1e-13)^2);
Efreq = zeros(1,Nt+1);
for i=1:(Nt+1)
    Efreq = runge_kuttaFT_E(Efreq,i);
end

%%
subplot(3,1,1)
plot(abs(p(:,Nt+1).^2))
head = sprintf('%s    |Pk|^2 at the end', date);
title(head)

subplot(3,1,2)
h = area(time, Et);
h.FaceColor = 0.8*[1,1,1];
h.LineStyle = 'None';
h.ShowBaseLine = 'off';
hold on

plot(time,abs(pt.^2)/max(abs(pt.^2)))

hold on

plot(time, sum(f)/max(sum(f)))
head = sprintf('%s    |P(t)|^2 ', date);
title(head)


subplot(3,1,3)
plot(time, Et)
head = sprintf('%s E(t)', date);
title(head)


figure
subplot(2,1,1)
plot(imag(pfreq))
title('Im[P({\omega} )]')

subplot(2,1,2)
plot(imag(pfreq)./abs(Efreq))
title( 'Im[P({\omega} )]/E({\omega})]')

