function [TIME, TE, NE, ND, ZEFF, SXR, R, Z, RAXIS, ZAXIS, ZEFF2D, NE2D, TE2D, SXR2D] = nc_sxr(filename, tp, plotyn);
% Function that simulates the SXR emission from
% TRANSP runs assuming that only braking radiation is
% contributing to the signal


% Example:
% filename = '/home/marco/Documents/MAST/TRANSP/RUNs/29976/U39/29976U39.CDF';
% tp = 0.2;
% plotyn = 1;
% [TE, NE, ND, ZEFF, SXR,  R, Z, RAXIS, ZAXIS, ZEFF2D, NE2D, TE2D, SXR2D] = nc_sxr(filename, tp, plotyn);


% Read the background plasma parameters
TIME = nc_read(filename, 'TIME');
TE = nc_read(filename, 'TE');
NE = nc_read(filename, 'NE');
ND = nc_read(filename, 'ND');
ZEFF = nc_read(filename, 'ZEFFP');
X = nc_read(filename, 'X');
XB = nc_read(filename, 'XB');
%[R, Z, RAXIS, ZAXIS] = nc_fluxsurfaces(filename);
[r, z, RAXIS, ZAXIS, PF, R, Z] = nc_fluxsurfaces(filename);

% Calculates the SXR signal
SXR = (NE.^2).*ZEFF.*sqrt(TE);  

% Create an empty array
[m, n] = size(SXR);
ZEFF_interp = zeros(m,n);
NE_interp = zeros(m,n);
TE_interp = zeros(m,n);
SXR_interp = zeros(m,n);

% Interpolates the data on the flux surface boundaries
for k = 1:m
  ZEFF_interp(k,:) = interp1(X(k,:), ZEFF(k,:), XB(k,:));
  NE_interp(k,:) = interp1(X(k,:), NE(k,:), XB(k,:));
  TE_interp(k,:) = interp1(X(k,:), TE(k,:), XB(k,:));
  SXR_interp(k,:) = interp1(X(k,:), SXR(k,:), XB(k,:));
end
u = isnan(ZEFF_interp);
ZEFF_interp(u) = 0;
clear u  
u = isnan(TE_interp);
TE_interp(u) = 0;
clear u  
u = isnan(NE_interp);
NE_interp(u) = 0;
clear u  
u = isnan(SXR_interp);
SXR_interp(u) = 0;
clear u  

n = n + 1;	% add one flux sruface for the magnetic axis
k = 50;
ZEFF2D = zeros(m, n, k);
NE2D = zeros(m, n, k);
TE2D = zeros(m, n, k);
SXR2D = zeros(m, n, k);

% the first flux surface contains the magnetic axis
for i = 1:m
  ZEFF2D(i,1,:) = ZEFF_interp(i,1);
  NE2D(i,1,:) = NE_interp(i,1);
  TE2D(i,1,:) = TE_interp(i,1);
  SXR2D(i,1,:) = SXR_interp(i,1);
end

% 2D profiles for all time intervals
for i = 1: m
	for j = 2:n
		ZEFF2D(i,j,:) = ZEFF_interp(i,j-1);
		NE2D(i,j,:) = NE_interp(i,j-1);
		TE2D(i,j,:) = TE_interp(i,j-1);		
		SXR2D(i,j,:) = SXR_interp(i,j-1);
	end
end
u = isnan(SXR2D);
SXR2D(u) = 0;
clear u  
u = isnan(ZEFF2D);
ZEFF2D(u) = 0;
clear u  
u = isnan(TE2D);
TE2D(u) = 0;
clear u  
u = isnan(NE2D);
NE2D(u) = 0;
clear u   
  

% plot the selected data

if (plotyn == 1)
  ip = max(find(TIME <= tp)); 
  % Read the flux surfaces 
  figure(1, 'position', [100 100 800 800])
    subplot(2,2,1)
      pcolor(squeeze(R(ip,:,:)), squeeze(Z(ip,:,:)), squeeze(ZEFF2D(ip,:,:))); 
      shading interp; 
      caxis([1, max(ZEFF(:))])
      colorbar
      hold on
      plot(RAXIS(ip) , ZAXIS(ip) , 'k+', 'linewidth', 2, 'markersize', 12)
      hold off
      axis equal
      xlabel('R (cm)')
      ylabel('Z (cm)')
      title(['Zeff @ ' num2str(TIME(ip)) ' s'])
     subplot(2,2,2)
      pcolor(squeeze(R(ip,:,:)), squeeze(Z(ip,:,:)), squeeze(NE2D(ip,:,:))); 
      shading interp;
      caxis([1, max(NE(:))])
      colorbar
      hold on
      plot(RAXIS(ip) , ZAXIS(ip) , 'k+', 'linewidth', 2, 'markersize', 12)
      hold off      
      axis equal
      xlabel('R (cm)')
      ylabel('Z (cm)')
      title(['Ne (cm-3) @ ' num2str(TIME(ip)) ' s']) 
    subplot(2,2,3)
      pcolor(squeeze(R(ip,:,:)), squeeze(Z(ip,:,:)), squeeze(TE2D(ip,:,:))); 
      shading interp; 
      caxis([0 max(TE(:))])
      colorbar
      hold on
      plot(RAXIS(ip) , ZAXIS(ip) , 'k+', 'linewidth', 2, 'markersize', 12)
      hold off      
      axis equal
      xlabel('R (cm)')
      ylabel('Z (cm)')
      title(['Te (eV) @ ' num2str(TIME(ip)) ' s'])
    subplot(2,2,4)
      pcolor(squeeze(R(ip,:,:)), squeeze(Z(ip,:,:)), squeeze(SXR2D(ip,:,:))); 
      shading interp; 
      caxis([0, max(SXR(:))])
      colorbar
      hold on
      plot(RAXIS(ip) , ZAXIS(ip) , 'k+', 'linewidth', 2, 'markersize', 12)
      hold off      
      axis equal
      xlabel('R (cm)')
      ylabel('Z (cm)')
      title(['SXR (a.u.) @ ' num2str(TIME(ip)) ' s'])     
endif
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



