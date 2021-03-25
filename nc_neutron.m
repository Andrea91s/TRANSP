function [TIME, TTNTX, BTNTX, BBNTX, THNTX, R, Z, X, XB, RAXIS, ZAXIS, NEUTT, MNEUT, BBNTS_DD, BTNTS_DD, NEUTX, DVOL] = nc_neutron(filename, plotyn)
         #[TIME, TTNTX, BTNTX, BBNTX, THNTX, R, Z, X, XB, RAXIS, ZAXIS, NEUTT, MNEUT, BBNTS_DD, BTNTS_DD, NEUTX, DV, UR, UZ, UTTNTX, UBTNTX, UBBNTX, UTHNTX, UDV]
% Script to access neutrons data from TRANSP output file
% INPUT
% filename	string		name of the NETCDF file containing TRANSP output
% plotyn 	integer		keyword: 0 no plot, 1 plot
%
% OUTPUT
% TIME		array		time in sec
% TTNTX		array 		total neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% BTNTX		array 		beam-thermal neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% BBNTX		array 		beam-beam neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% THNTX		array 		thermal neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% R		array 		flux surface radial coordinate in cm
% Z		array 		flux surface vertical coordinate in cm
% X		array		flux surfaces centres
% XB		array		flux surfaces boundaries
% RAXIS		array 		magnetic axis radial coordinate in cm
% ZAXIS 	array 		magnetic axis vertical coordinate in cm
% NEUTT		array 		total neutron yield in s-1
% MNEUT		array 		measure total neutron yield
% BBNTS_DD	array 		DD beam-beam neutron yield s-1
% BTNTS_DD	array 		DD beam-thermal neutron yield s-1
% NEUTX		array 		DD thermal neutron yield in s-1
% 
% Example:
%filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/K25/C12/99999C12.CDF';
%filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/K26/99999K26.CDF';
% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29881/29881O07_fi_1.cdf'
% nc_neutron(filename, 1)


% Read the time
TIME = nc_read(filename, 'TIME');   % time in sec
NTHETA=50;
method=1;

% Read the total neutron emissivities defined on the centres of
% the flux surfaces boundaries
TTNTX = nc_neutronemissivity(filename, 'TTNTX', NTHETA, method);
BTNTX = nc_neutronemissivity(filename, 'BTNTX', NTHETA, method);
BBNTX = nc_neutronemissivity(filename, 'BBNTX', NTHETA, method);
THNTX = nc_neutronemissivity(filename, 'THNTX', NTHETA, method);
DVOL = nc_read(filename, 'DVOL');

% Read the flux surfaces boundaries, the normalized boundaries and centres
[R, Z] = nc_fluxsurfaces(filename);
X = nc_read(filename, 'X');
XB = nc_read(filename, 'XB');

% Read the magnetic axis position
RAXIS = nc_read(filename, 'RAXIS');
ZAXIS = nc_read(filename, 'YAXIS');

% Read the total neutron yields
%YNM = ncread(filename, 'MNEUT');   % measured neutrons/s
NEUTT = nc_read(filename, 'NEUTT');   % TRANSP calculated neutrons/s

% Read the neutron yield components
BBNTS_DD = nc_read(filename, 'BBNTS_DD');	% DD beam-beam neutrons/s
BTNTS_DD = nc_read(filename, 'BTNTS_DD');	% DD beam-target neutrons/s
NEUTX = nc_read(filename, 'NEUTX');		% thermal neutrons/s

% Calculates the total neutron yield from the emissivities



% Plot the time traces
if (plotyn == 1)
  figure(1)
    plot(TIME, NEUTT, 'r', 'linewidth', 2, TIME, BBNTS_DD, 'b', 'linewidth', 2, TIME, BTNTS_DD, 'm', 'linewidth', 2, TIME, NEUTX, 'g', 'linewidth', 2) 
    axis([TIME(1) TIME(end) 0 1.05*max(NEUTT)])
    xlabel('time (s)')
    ylabel('Neutron Yield (1/s)')
    legend('Total', 'BB', 'BT', 'TH', 'location', 'northwest')
    title(strrep(filename, '_', ' '));
endif

close all;

time = 7.54;
idx = max(find(TIME <= time));

TH = squeeze(THNTX(idx,:,:));
BT = squeeze(BTNTX(idx,:,:));
BB = squeeze(BBNTX(idx,:,:));
r = squeeze(R(idx,:,:));
z = squeeze(Z(idx,:,:));

figure(1)
subplot(2,2,1)
pcolor(r,z,TH)
shading interp
colorbar
caxis([0 8e8])
title('thermal')

subplot(2,2,2)
pcolor(r,z,BT)
shading interp
colorbar
caxis([0 8e8])
title('beam-thermal')

subplot(2,2,3)
pcolor(r,z,BB)
shading interp
colorbar
caxis([0 8e8])
title('beam-beam')

subplot(2,2,4)
pcolor(r,z,TH+BT+BB)
shading interp
colorbar
title('total')










	
