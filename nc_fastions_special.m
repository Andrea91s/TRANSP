function [R, Z, FIDD] = nc_fastions_special(filename, plotyn, saveyn)
% Script to access fast ions data from the special TRANSP output file
%
% INPUT
% filename = '/home/andrea/Documents/TRANSP/29880/29880U16_fi_3.cdf';
% filename	string		name of the special TRANSP output
% plotyn	integer		0 no plot, 1 plot
% saveyn	integer		0 no save, 1 save
%
% OUTPUT
% R		array		radial coordinates in a given time in cm
% Z		array		vertical coordinates in a given time in cm
% FIDD		array		energy/pitch angle integrated FID



% Read the energy and pitch angles
PA = nc_read(filename, 'A_D_NBI');
EN = nc_read(filename, 'E_D_NBI');

% Read the coordinates
R2D = nc_read(filename, 'R2D');
Z2D = nc_read(filename, 'Z2D');

% Read the fast ion distribution function
FID = nc_read(filename, 'F_D_NBI');

% Integrate the FI over energy and pitch angle
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);
w = sum(sum(FID,3),2)*DPA*DEN;

% Generate the 2D map
[R Z] = meshgrid (linspace(0,200,200), linspace(-100, 100, 200));
FIDD = griddata(R2D, Z2D, w, R, Z, 'linear');
u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u

% Plot the FID distribution
if (plotyn == 1)
  %load mycolormaps     
    figure(2)
	pcolor(R, Z, FIDD);
	shading interp
	%set(gcf,'Colormap',MASTcmap) 
	xlabel('R (cm)')
	ylabel('Z (cm)')
	axis([0 200 -100 100])
	axis equal
	title(['Fast Ion Distribution'])
	colorbar
  
endif

% Save the data
if (saveyn == 1)
  save(['TRANSP_FID.dat'], 'FIDD', '-ascii')
endif