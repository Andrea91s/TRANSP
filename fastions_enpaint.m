function [FIDS, FIDn, FIDN, FIDT, FIDE, FIDP, id] = fastions_enpaint(R2D, Z2D, PA, EN, FID, BMVOL, EI1, EI2, PI1, PI2, plotyn, RR, RZ, CM, IE)
% Script that calculates the fast ions density in a point on the poloidal cross-section
% and in a user specifed energy region and range
% INPUT
% R2D		array		radial position at which the FID has been calculated in cm
% Z2D		array		vertical position at which the FID has been calculated in cm
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% BMVOL		array		zone volumes
% EI1		float		lower integral limit in energy
% EI2		float		upper integral limit in energy
% PI1		float		lower integral limit in pitch angle
% PI2		float		upper integral limit in pith angle
% plotyn	integer		0 = no plot, 1 = polt
% RR		real		Requested radial coordinate to plot the FID
% RZ		real		Requested vertical position to plot the FID	
% CM 		real		Z max of the FID: 0 for auto
% IE		real		maximum injection energy for plotting
%
% OUTPUT
% FIDS		array		fast ion density in the energy/pitch angle range
% FIDn		array		local (R,Z) fast ion density (cm-3) integrated over the specified region 
% FIDN		float		local (R,Z) total number of fast ions integrated over th specified region
% FIDT		float		total number of fast ions integrated over th specified region
% FIDE		array		fast ion distribution function integrated over all pitch angles as a function of energy
% FIDP		array		fast ion distribution function integrated over all energy as a function of pitch ange
% id		integer		array index of FID which is closest to the requested point in the poloidal plane
% Example
% fastions_enpaint(R2D, Z2D, PA, EN, FID, BMVOL, EI1, EI2, PI1, PI2, 0, 90, 0, 0, 0);



% Create the array for plotting the 2D FID
[pa, en] = meshgrid(PA, EN);
[X Y] = meshgrid (linspace(0,200,200), linspace(-100, 100, 200));


% Check the integration limits
if (EI1 < min(EN))
  EI1 = min(EN);
endif
if (PI1 < min(PA))
  PI1 = min(PA);
endif


% Find the indexes of the integral limits in the array
i1 = max(find(EN <= EI1))
  ENS1 = EN(i1);
i2 = max(find(EN <= EI2));
  ENS2 = EN(i2);
i3 = max(find(PA <= PI1))+1
  PAS1 = PA(i3);
i4 = max(find(PA <= PI2));
  PAS2 = PA(i4);

  
% Calculates the fast ion density (cm-3) and the fast ion numbers for all (R,Z) 
% in the specified energy/pitch angle range
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);
FIDn = sum(sum(FID(:,i3:i4,i1:i2),3),2)*DPA*DEN;
FIDN = sum((sum(FID(:,i3:i4,i1:i2),2)*DPA)*DEN,3).*BMVOL;
FIDT = sum(sum((sum(FID(:,i3:i4,i1:i2),2)*DPA)*DEN,3).*BMVOL,1);

% Interpolates the fast ion density (cm-3) FIDn over the X,Y cooridnates
% to make a plot of the FID density integrated in the energy and pitch angle
% selected regions as a function of (R,Z)
FIDni = griddata(R2D, Z2D, FIDn, X, Y, 'linear');
u = find(isnan(FIDni) == 1);
FIDni(u) = 0;
clear u

% Caclualtes the fast ion distribution function limited to the specified energ/pitch angle range
% for all (R,Z)
FIDS = zeros(size(FID));
FIDS(:,i3:i4, i1:i2) = FID(:,i3:i4, i1:i2);

% Calculates the FID(E) integrated over the pitch angle and FID(PA) integrated over the energy
FIDE = squeeze(sum(FIDS,2))*DPA;
FIDP = squeeze(sum(FIDS,3))*DEN;

% Output the results
printf('Total number of fast ions in the range %f to %f keV and %f to %f\n', EI1/1000, EI2/1000, PI1, PI2);
printf('From FBM: %1.4g\n', FIDT);


% Find the index of the FID which is closest in space
% to the required point
if ((RR == 0) && (RZ == 0))
  id = 1;
else
  d = sqrt((R2D-RR).^2 + (Z2D-RZ).^2);
  di = find(d == min(d));
  id = di(1);
endif
if (CM == 0)
  CM = max(max(FID(id,:,:)));
endif


% Plot the time traces
if (plotyn == 1)
  figure(5, 'position', [100 100 800 800])
  subplot(2,2,1)
  plot(R2D, Z2D, '.', R2D(id), Z2D(id), 'ro', RR, RZ, 'ks')
   xlabel('R (cm)')
   ylabel('Z (cm)')
   legend('FID points', 'Closest FID', 'Requested FID')
   axis([0 200 -100 100])
   axis equal
   
  subplot(2,2,2)
    pcolor(pa, en, squeeze(FIDS(id,:,:))');
    shading interp
    xlabel('pitch angle')
    ylabel('energy eV')
    axis([-1 1 0 IE 0 CM])
    caxis([0 CM])
    %axis tight
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV, pitch angle ' num2str(PAS1) ' - ' num2str(PAS2) '.'])
    colorbar
 
 subplot(2,2,3) 
    plot(EN/1000, FIDE(id,:), 'linewidth', 2)
    xlabel ('Energy (keV)')
    ylabel ('dFID(E)/dE density')
    axis([0 IE/1000 0 1.1*max(FIDE(id,:))])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])    
  
  subplot(2,2,4) 
    plot(PA, FIDP(id,:), 'linewidth', 2)
    xlabel ('pitch angle')
    ylabel ('dFID(Pa)/dPa density')
    axis([-1 1 0 1.1*max(FIDP(id,:))])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])    

  figure(6)  
  load mycolormaps     
  pcolor(X, Y, FIDni);
    shading interp
    set(gcf,'Colormap',MASTcmap) 
    xlabel('R (cm)')
    ylabel('Z (cm)')
    axis([40 160 -70 70])
    axis equal
    title(['Fast Ion Distribution - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV, pitch angle ' num2str(PAS1) ' - ' num2str(PAS2) '.' ])
    colorbar    
endif














	
