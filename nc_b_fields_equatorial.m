function [RR, TIME, BP, BT, BTOT] = nc_b_fields_equatorial(filename)
% Script to read the magnetic fields components (BT and BP) along the 
% major radius, that is for Z = 0
%
% INPUT
% filename		name of the TRANSP file containing the results
%
% OUTPUT
% R		radius of the flux surfaces
% TIME		time array
% BT		toroidal B field
% BP		poloidal B field
% BTOT		total B field

% Read the arrays
RMAJM = nc_read(filename, 'RMAJM');
RAXIS = nc_read(filename, 'RAXIS');
BTX = nc_read(filename, 'BTX');
TIME = nc_read(filename, 'TIME');
FBTX = nc_read(filename, 'FBTX');
FBPBT = nc_read(filename, 'FBPBT');


% Calculates BT
BT = BTX.*FBTX;

% Calculates BP
BP = BT.*FBPBT;

% Calculates Btotal
BTOT = sqrt(BT.^2 + BP.^2);

% Calculates the total field
BTOT = sqrt(BT.^2 + BP.^2);


% Calculates the actual radial position where
% the magnetic fields are calculates
RR = zeros(size(RMAJM));
for k = 1:length(RAXIS)
  RR(k,:) = RAXIS(k) + RMAJM(k,:) - RMAJM(k,61);
end


