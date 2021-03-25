% Routine to test if I use correctly the neutron emissivity
% 
% The idea is to read NEUTT and compare it with sum of the
% neutron emissivity that I calculate and then use to
% estimate the expected count rates at the detectors of
% the neutron camera.


% Specify the file to read:
filename = '/home/andrea/Desktop/D02.CDF';

% Read the data in the usual way (but with method = 1 that is without interpolating the point
% on the magnetic axis):
[TIME, TTNTX, BTNTX, BBNTX, THNTX, R, Z, X, XB, RAXIS, ZAXIS, NEUTT, MNEUT, BBNTS_DD, BTNTS_DD, NEUTX, DVOL] = nc_neutron(filename, plotyn = 1);

% Number of elements
[i,j,k] = size(TTNTX);


% Calculate the sum of the thermal, beam-thermal and beam-beam components from the
% neutron emissivity;
W = BBNTX + BTNTX + THNTX;

% Since some of the values are NAN I set them to zero (these come form the linear interpolation
% on BX boundaries:
W(find(isnan(W))) = 0;
TTNTX(find(isnan(TTNTX))) = 0;

% Calculates the total neutron yield as the sum of the three components ...
YN = W(:,1:j,1).*DV;
YN = sum(sum(YN,3),2);

% ... and from the total value read from the TRANSP output:
TYN = TTNTX(:,1:j,1).*DV;
TYN = sum(sum(TYN,3),2);

% ... and the one used to make plots with pcolor
UT = UTTNTX (1:i,1:j+1,1).*UDV(1:i,1:j+1);
UT = sum(sum(UT,3),2);

% Compare YN with NEUTT
figure(1)
plot(TIME,NEUTT, 'linewidth', 2, TIME, YN, 'r--', 'linewidth', 2, TIME, TYN, 'k+', 'linewidth', 2, TIME, UT, 'om', 'linewidth', 2)
xlabel('time (s)', 'fontsize', 14)
ylabel('Neutron Yield (s-1)', 'fontsize', 14)
title(filename)
legend('NEUTT', 'YN', 'TYN', 'UT')

return
% Plots the emissivity profiles for a given time
nc_neutronemissivity_plot(TIME, UR, UZ, UTTNTX, UDV, time = 0.234, fs = 1);
nc_neutronemissivity_onaxis(TIME, R, Z, TTNTX, RAXIS, ZAXIS, time = 0.234, plotyn = 1);



% Calculates the neutron emissivity at a selected time (the first sawtooth data from the special output):
[r, z, ett, ebt, ebb, eth, RMAV, ZMAV, ri, zi, etti, NEUTTAV, MNEUTAV, DVAV, YNT] = nd_time_average(ts = 0.233, te = 0.235, TIME, R, Z, TTNTX, BTNTX, BBNTX, THNTX, DV, NEUTT, MNEUT, RAXIS, ZAXIS, plotyn = 1, saveyn = 0, '');


% Read the special output

filename_neut = '/home/andrea/Documents/TRANSP_FID/29881/29881O07_neut_1.cdf';
filename_fi = '/home/andrea/Documents/TRANSP_FID/29881/29881O07_fi_1.cdf';
[RNA, ZNA, BBN4, BTN4, BMVOL, YNTA] = nc_neutron_special(filename_neut, filename_fi);

% Combines the special output and the regular outputs
[ett_na, ebt_na, ebb_na, etti_na, YNTA, r_zero, ett_na_zero] = nd_special_emissivity(r, z, eth, DVAV, RNA, ri, zi, ZNA, BBN4, BTN4, BMVOL, plotyn = 1, saveyn = 0, '');

% Evaluates the asymmetry
[e_asym, ei_asym] = nd_emissivity_asymmetry(r, z,  ett, ett_na, ri, zi, etti, etti_na, plotyn = 1, saveyn = 0);
