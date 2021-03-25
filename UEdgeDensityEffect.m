filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U01/29904U01.CDF';
T = nc_read(filename, 'TIME');
U01_Yn = nc_read(filename, 'NEUTT');

filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904//U06/29904U06.CDF';
U06_Yn = nc_read(filename, 'NEUTT');

plot(T, U01_Yn, 'k', 'linewidth', 2, T, U06_Yn, 'r', 'linewidth', 2)
xlabel('time (s)', 'fontsize', 12)
ylabel('NEUTT (s^{-1})', 'fontsize', 12)
title('Plasma Discharge 29904', 'fontsize', 14)
legend('Run U01', 'Run U06', 'location', 'northwest')


filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U01/F29904.RCY';
UFILE01 = ufile_read_neutral(filename ,0);

filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U06/F29904.RCY';
UFILE06 = ufile_read_neutral(filename ,0);

figure
semilogy(UFILE01.t, UFILE01.data, 'k', 'linewidth', 2, UFILE06.t, UFILE06.data, 'r', 'linewidth', 2)
xlabel('time (s)', 'fontsize', 12)
ylabel('Edge Neutral Flux (s^{-1})', 'fontsize', 12)
title('Plasma Discharge 29904 - Edge Recycling', 'fontsize', 14)
legend('Run U01', 'Run U06', 'location', 'northwest')


filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U01/29904U01.CDF';
DN0WD = nc_read(filename, 'DN0WD');
filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29904/U06/29904U06.CDF';
DN0WDU06 = nc_read(filename, 'DN0WD');
DN0WDU01 = DN0WD;

figure(1, 'position', [100 100 1400 500])
subplot(1,2,1)
pcolor(log10(DN0WDU01)); 
shading flat; 
colorbar
caxis([min(log10(DN0WDU01(:))) max(log10(DN0WDU06(:)))])
xlabel('r/a (a.u.)', 'fontsize', 12)
ylabel('time (a.u.)', 'fontsize', 12)
title ('U01 - DN0WD (cm^{-3})', 'fontsize', 14)

subplot(1,2,2)
pcolor(log10(DN0WDU06)); shading flat; 
colorbar
caxis([min(log10(DN0WDU01(:))) max(log10(DN0WDU06(:)))])
xlabel('r/a (a.u.)', 'fontsize', 12)
ylabel('time (a.u.)', 'fontsize', 12)
title ('U06 - DN0WD (cm^{-3})', 'fontsize', 14)

figure(2)
hold on
h = semilogy(1:60, [DN0WDU01(100,:)' DN0WDU06(100,:)'], {'k', 'r'});
set (h, 'linewidth', 2)
h = semilogy(1:60, [DN0WDU01(70,:)' DN0WDU06(70,:)'], {'ko', 'g'});
set (h, 'linewidth', 2)
hold off
xlabel('r/a (a.u.)', 'fontsize', 12)
ylabel('DN0WD (cm^{-3})', 'fontsize', 12)
legend('U01: t(100)', 'U06: t(100)', 'U01:t(70)', 'U06:t(70)')

return

semilogy(1:60, [DN0WDU01(70,:)' DN0WDU06(70,:)'])
semilogy(1:60, [DN0WDU01(100,:) DN0WDU06(100.:)])
semilogy(1:60, [DN0WDU01(100,:) DN0WDU06(100.:)])
semilogy(1:60, [DN0WDU01(100,:) DN0WDU06(100,:)])
semilogy(1:60, [DN0WDU01(100,:) DN0WDU06(100,:)])
semilogy(1:60, [DN0WDU01(100,:)' DN0WDU06(100,:)'])
semilogy(1:60, [DN0WDU01(75,:)' DN0WDU06(75,:)'])
semilogy(1:60, [DN0WDU01(70,:)' DN0WDU06(70,:)'])

