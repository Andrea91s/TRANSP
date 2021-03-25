% Routine to test if I use correctly the neutron emissivity
% 
% The idea is to read NEUTT and compare it with sum of the
% neutron emissivity that I calculate and then use to
% estimate the expected count rates at the detectors of
% the neutron camera.

% Input data
% %{
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '29904';
transp_id = '/29904O04';
filename = [directory pulse transp_id '.CDF'];
time = [0.194 0.213 0.221 0.241 0.256];
ts = [0.193 0.212 0.220 0.240 0.255];
te = [0.195 0.215 0.222 0.243 0.257];
%}

%{
% Series 29210
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '29210';
transp_id = '/O22/29210O22';
filename = [directory pulse transp_id '.CDF'];
te = [0.135,0.165,0.195,0.220,0.240,0.255,0.265];
ts = te - 0.003;
time = te - 0.0015;
%}


%%{
% Series 29924
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '29924';
transp_id = '/O18/29924O18';
filename = [directory pulse transp_id '.CDF'];
te = [0.185, 0.205, 0.235, 0.265, 0.295, 0.305, 0.315, 0.325];
ts = te - 0.003; 
time = te - 0.0015;
%}


%{
% Series 27938
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '27938';
transp_id = '/U11/27938U11';
filename = [directory pulse transp_id '.CDF'];
te = [0.252, 0.27065, 0.282, 0.315, 0.340, 0.369, 0.378];
ts = te - 0.005; 
time = te - 0.0015;
%}



%{
% Series 29976
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '29976';
transp_id = '/O69/29976O69';
filename = [directory pulse transp_id '.CDF'];
te = [0.180, 0.201, 0.206, 0.214, 0.219, 0.226, 0.233];
ts = te - 0.003; 
time = te - 0.0015;
%}

%{
% Series 29880
directory = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/';
pulse = '29880/U04/29880';
transp_id = 'U04';
filename = [directory pulse transp_id '.CDF'];
%te = [0.235 0.258 0.279 0.295];
%ts = [0.232 0.255 0.276 0.292];
ts = [0.239 0.260 0.282 0.300];
te = [0.242 0.263 0.285 0.303];
time = te - 0.0015;

selected_idx = 1;
fi_id = [transp_id '_fi_' num2str(selected_idx) '.cdf'];
neut_id = [transp_id '_neut_' num2str(selected_idx) '.cdf'];
filename_neut = [directory pulse neut_id];
filename_fi =  [directory pulse fi_id];

time = time(1);
ts = ts(1);
te = te(1);
%}


% Selected data
selected_idx = 3;
fi_id = [transp_id '_fi_' num2str(selected_idx) '.cdf'];
neut_id = [transp_id '_neut_' num2str(selected_idx) '.cdf'];
filename_neut = [directory pulse neut_id];
filename_fi =  [directory pulse fi_id];

time = time(selected_idx);
ts = ts(selected_idx);
te = te(selected_idx);



% Read the data in the usual way (but with method = 1 that is without interpolating the point
% on the magnetic axis):
[TIME, TTNTX, BTNTX, BBNTX, THNTX, R, Z, X, XB, RAXIS, ZAXIS, NEUTT, MNEUT, BBNTS_DD, BTNTS_DD, NEUTX, DV, UR, UZ, UTTNTX, UBTNTX, UBBNTX, UTHNTX, UDV] = nc_neutron(filename, plotyn = 0);

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

plotyn = 1;
% Compare YN with NEUTT
if (plotyn == 1)
figure(1)
plot(TIME,NEUTT, 'linewidth', 2, TIME, YN, 'r--', 'linewidth', 2, TIME, TYN, 'k+', 'linewidth', 2, TIME, UT, 'om', 'linewidth', 2)
xlabel('time (s)', 'fontsize', 14)
ylabel('Neutron Yield (s-1)', 'fontsize', 14)
title(filename)
legend('NEUTT', 'YN', 'TYN', 'UT')
endif

% Plots the emissivity profiles for a given time
if (plotyn == 1)
nc_neutronemissivity_plot(TIME, UR, UZ, UTTNTX, UDV, time, fs = 1);
nc_neutronemissivity_onaxis(TIME, R, Z, TTNTX, RAXIS, ZAXIS, time, plotyn = 1);
endif


% Calculates the neutron emissivity at a selected time (the first sawtooth data from the special output):
[r, z, ett, ebt, ebb, eth, RMAV, ZMAV, ri, zi, etti, NEUTTAV, MNEUTAV, DVAV, YNT] = nd_time_average(ts, te, TIME, R, Z, TTNTX, BTNTX, BBNTX, THNTX, DV, NEUTT, MNEUT, RAXIS, ZAXIS, plotyn, saveyn = 0, '');


% Read the special output

[RNA, ZNA, BBN4, BTN4, BMVOL, YNTA] = nc_neutron_special(filename_neut, filename_fi);

% Combines the special output and the regular outputs
[ett_na, ebt_na, ebb_na, etti_na, YNTA, r_zero, ett_na_zero] = nd_special_emissivity(r, z, eth, DVAV, RNA, ri, zi, ZNA, BBN4, BTN4, BMVOL, plotyn, saveyn = 0, '');

% Evaluates the asymmetry
[e_asym, ei_asym] = nd_emissivity_asymmetry(r, z,  ett, ett_na, ri, zi, etti, etti_na, plotyn, saveyn = 0);
