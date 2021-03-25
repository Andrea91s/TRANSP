function [R, Z, RAXIS, ZAXIS, PFLUX, UR, UZ] = nc_fluxsurfaces(filename);
% Function that returns the flux surfaces from the NetCDF file
% specified in the
%
% INPUT
% filename		NetCDF file name
%
% OUTPUT
% R				time x R x theta coordinates in cm
% Z				time x Z x theta coordinates in cm
% RAXIS
% ZAXIS
% PFLUX         poloidal flux 
% UR            time x R x theta coordinates in cm with added on-axis fake surface
% UZ            time x Z x theta coordinates in cm with added on-axis fake surface


% Number of poloidal angles: note that when this script is called in conjuction with
%                            with the script nc_neutronemissivity.m the number of
%                            poloidal angles must be set to the same value on both
%                            scripts!
NTHETA = 360;
NTHETA = 50;

% Read the moments
RMC00 = nc_read(filename, 'RMC00');
YMC00 = nc_read(filename, 'YMC00');

% Read the poloidal flux on the boundary XB
PF = nc_read(filename, 'PLFLX');

% Create the arrays
[M, N] = size(RMC00);
P = 10;
RMC = zeros(M,N,P);
RMS = zeros(M,N,P);
YMC = zeros(M,N,P);
YMS = zeros(M,N,P);

% Read the data into the arrays
for p = 1:P
  if (p < 10)
    SRC = ['RMC0' num2str(p)];
    SRS = ['RMS0' num2str(p)];
    SYC = ['YMC0' num2str(p)];
    SYS = ['YMS0' num2str(p)];
  else
    SRC = ['RMC' num2str(p)];
    SRS = ['RMS' num2str(p)];
    SYC = ['YMC' num2str(p)];
    SYS = ['YMS' num2str(p)];   
  endif
  RMC(:,:,p) = nc_read(filename, SRC);
  RMS(:,:,p) = nc_read(filename, SRS);
  YMC(:,:,p) = nc_read(filename, SYC);
  YMS(:,:,p) = nc_read(filename, SYS);
end

% Read the magnetic axis position
RAXIS = nc_read(filename, 'RAXIS');
ZAXIS = nc_read(filename, 'YAXIS');

% Poloidal angle
theta = linspace(0,2*pi, NTHETA);
K = length(theta);

% Create the array for the coordinates of
% the flux surfaces
R = zeros(M, N, K);
Z = zeros(M, N, K);
PFLUX = zeros(M, N, K); 

% Define the cosine and sine moments
CM = zeros(P, K);
SM = zeros(P, K);
for p = 1:P
  CM(p,:) = cos(p*theta);
  SM(p,:) = sin(p*theta);
end

% Multiplies the moments and the cosines, sines
for m = 1:M
  R(m,:,:) = squeeze(RMC(m,:,:))*CM + squeeze(RMS(m,:,:))*SM;
  Z(m,:,:) = squeeze(YMC(m,:,:))*CM + squeeze(YMS(m,:,:))*SM;
end

% Generates the coordinates without the fake flux surface at
% the magnetic axis
for n = 1:N
  for k = 1:K
    R(:,n,k) = RMC00(:,n) + squeeze(R(:,n,k)); 
    Z(:,n,k) = YMC00(:,n) + squeeze(Z(:,n,k)); 
  end
end

% Generates the poloidal flux
for n = 1:N
    PFLUX(:,n,:) = PF(:,n)*ones(1,NTHETA);
end

% Generates the coordinates adding a fake  flux surface
% on the magnetic axis

UR = zeros(M, N+1, K);
UZ = zeros(M, N+1, K);

for k = 1:K
  UR(:,1,k) = RAXIS(:); 
  UZ(:,1,k) = ZAXIS(:); 
end

UR(:,2:end,:) = R;
UZ(:,2:end,:) = Z;

%{
clear R Z
R = UR;
Z = UZ;
clear UR UZ
%}  
  
%%%%%%%%%%%%%%%%%%%% obtain the normalized poloidal flux

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



