
addpath('/home/andrea/Documents/MASTOrbit');
[EFIT] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/g029909.00216',1);
[TEQ] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/29909EQDSK.eqdsk',1); #THIS IS FROM TRXPL
close all;
figure(1)
contour(EFIT.R, EFIT.Z, EFIT.PSI_2D, 50,'k', 'linewidth',2)
hold all;
contour(TEQ.R, TEQ.Z, TEQ.PSI_2D, 50,'r--','linewidth',2)
plot(EFIT.R_axis, EFIT.Z_axis,'+k', 'markersize', 20, 'linewidth', 3)
plot(TEQ.R_axis, TEQ.Z_axis, '+r', 'markersize', 20, 'linewidth', 3)
axis equal
xlim([min(TEQ.R), max(TEQ.R)])
ylim([min(TEQ.Z), max(TEQ.Z)])
legend('EFIT', 'TEQ')
set (gca, 'fontsize', 14)
xlabel('R (m)','fontsize', 14)
ylabel('Z (m)','fontsize', 14)
title('Poloidal Flux \Psi_{p}', 'fontsize', 14, 'fontweight', 'bold')

close all;

figure(2)
  % Plasma pressure
subplot(3,1,1)
plot(EFIT.psi_rho_pol, EFIT.pressure, 'k', 'linewidth', 2);
hold all;
plot(TEQ.psi_rho_pol, TEQ.pressure, 'r', 'linewidth', 2);
xlabel('\rho_{pol}','fontsize', 14)
ylabel('Pressure (Newton/m^2)','fontsize', 14)
legend('EFIT', 'TEQ')
set (gca, 'fontsize', 14)
    
  % Safety factor
subplot(3,1,2)
plot(EFIT.psi_rho_pol, EFIT.q, 'k', 'linewidth', 2);
hold all;
plot(TEQ.psi_rho_pol, TEQ.q, 'r', 'linewidth', 2);
xlabel('\rho_{pol}','fontsize', 14)
ylabel('Safety factor','fontsize', 14)
legend('EFIT', 'TEQ')
set (gca, 'fontsize', 14)
    
  % T factor
subplot(3,1,3)
plot(EFIT.psi_rho_pol, EFIT.F, 'k', 'linewidth', 2);
hold all;
plot(TEQ.psi_rho_pol, TEQ.F, 'r', 'linewidth', 2);
xlabel('\rho_{pol}','fontsize', 14)
ylabel('Poloidal current function F (m Tesla)','fontsize', 14)      
legend('EFIT', 'TEQ')
set (gca, 'fontsize', 14)

close all;

figure (3)

    subplot(2,1,1)
    plot(EFIT.R, EFIT.psi_1D, 'k', 'linewidth', 2);
    hold all;
    plot(TEQ.R, TEQ.psi_1D, 'r', 'linewidth', 2);
    xlabel('R (m)','fontsize', 14)
    ylabel('stream function \Psi (Weber/rad)','fontsize', 14)     
    xlim([min(TEQ.R), max(TEQ.R)])
    set (gca, 'fontsize', 14)

    
    subplot(2,1,2)
        plot(EFIT.R, EFIT.B.POL.R_radial, 'k', 'linewidth', 2, ...
         EFIT.R, EFIT.B.POL.Z_radial, 'k--', 'linewidth', 2, ...
         EFIT.R, EFIT.B.TOR_radial, 'k-.', 'linewidth', 2)
         hold all;
         plot(TEQ.R, TEQ.B.POL.R_radial, 'r', 'linewidth', 2, ...
         TEQ.R, TEQ.B.POL.Z_radial, 'r--', 'linewidth', 2, ...
         TEQ.R, TEQ.B.TOR_radial, 'r-.', 'linewidth', 2)
    hl = legend('EFIT, B_{\theta,R}(R,Z_{0})', 'EFIT, B_{\theta,Z}(R,Z_{0})', 'EFIT, B_{\phi}(R,Z_{0})',
                'TEQ, B_{\theta,R}(R,Z_{0})', 'TEQ, B_{\theta,Z}(R,Z_{0})', 'TEQ, B_{\phi}(R,Z_{0})');
    set(hl, 'box', 'off', 'location', 'southwest', 'fontsize', 10)
    xlabel('R (m)','fontsize', 14)
    ylabel('Magnetic Field (T)','fontsize', 14)
    xlim([-0.1, max(TEQ.R)])
    ylim([-2, 2.5])    
    set (gca, 'fontsize', 14)

