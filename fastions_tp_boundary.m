function [X, BMIN, BMAX] = fastions_tp_boundary(filename, tb, R2D, RR, plotyn, saveyn)
% Script to access neutrons data from TRANSP output file
% INPUT
% filename	string		name of the NETCDF file containing TRANSP output
% tb		float		time of the required boundary
% R2D		array		radial coordinates from the fast ion distribution function
% RR		array		requested RR coordinate
%
% OUTPUT
% X		array		flux surface normalized coordinate
% BMIN		array		min B field on the flux surface
% BMAX		array		max B field on the flux surface
% TPB		array		passing - trapped boundary 
%
% Example
% filename = '/home/andrea/Documents/TRANSP/29880/29880U16.cdf';
% [R2D, Z2D, PA, EN, FID, X, Y, FD, FIDD, r, er] = fastions(filename, tb ,R2D, RR, 1, 0);


% Read the pitch angle scattering values
%TIME = nc_read(filename, 'TIME3');
BMIN = nc_read(filename, 'BMIN');
BMAX = nc_read(filename, 'BMAX');
X = nc_read(filename, 'X');


% Evaluate the passing - trapped boundary
TPB = sqrt(1-BMIN./BMAX);


% Find the required time
it = max(find(TIME <= tb));


% Find the required normalized flux surface
RM = mean(R2D(1:4)); 		% radial coordinate of the magnetic axis
if (RR >= RM)
  RE = abs(max(R2D) - RM);		% plasma minor a
else
  RE = abs(min(R2D) - RM);
end
xr = abs(RR-RM)/RE;


% Find the normalized flux surface corresponding to the requested one
ir = max(find(X(it,:) <= xr));
printf('The TP boundary at t = %f for r/a = %f is %f\n', TIME(it), X(it,ir), TPB(it,ir));


% Plot the boundary
if (plotyn == 1)
  TT = zeros(size(TPB));
  for k = 1:60
    TT(:,k) = TIME;
  end
  
  pcolor(TT, X, TPB)
  shading interp
  caxis([0 1])
  colorbar
  caxis([0 1])
  axis tight
  xlabel('time (s)')
  ylabel('Flux surface')
endif











	
