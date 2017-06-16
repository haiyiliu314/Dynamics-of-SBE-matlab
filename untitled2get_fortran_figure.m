% load('/home/liuhai/matlab_code/dynamics_coulomb/result/Run_Ben.mat')
% uiopen('/home/liuhai/matlab_code/dynamics_coulomb/Efreq.dat',1)
% freqgrid =(E001'*Ebind + 4*bind)/hbar;  
clear pt
pfreq1 = zeros(1, Nfreq+1);
test = importdata('/home/liuhai/dynamics/with_coulomb_Fortran/pt.dat', ',');
pt = (test(:,1)+1i*test(:,2))';
for i=1:(Nt)
    pfreq1 = runge_kuttaFT(pfreq1,i);
end
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
plot(freqgrid*hbar, imag(pfreq1))
title('Im[P({\omega})]')

subplot(5,1,2)
plot(freqgrid*hbar, imag(pfreq1./(Efreq))/max(abs(imag(pfreq1./(Efreq)))))
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


figure
plot(E001'*Ebind+4*Ebind, abs((imag(pfreq1./(Efreq))/max(abs(imag(pfreq1./(Efreq)))) - E1'/max(E1))./(E1'/max(E1)))*100)
xlabel('h\omega -(E_g+E_{1s})[meV]')
head = sprintf('%s Comparison on susceptibility', date);
title(head)