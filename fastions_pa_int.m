function [FIDDST, FIDDSP, FIDDSTT, FIDDSPT, id] = fastions_pa_int(R2D, Z2D, PA, EN, FID, EI1, EI2, filename, tr, RR, RZ, plotyn, CM, IE)
% Script that calculates for each point in the poloidal projection the trapped-passing boundary and then calculates
% the density of passing and trapped particles in all points in a given energy range.
% 
% NB! 	The trapped/passing boundary is a very crude approx. based on the assumption 
%	that the magnetic moment is conserved which is not in ST!! 
%
% INPUT
% R2D		array		radial position at which the FID has been calculated in cm
% Z2D		array		vertical position at which the FID has been calculated in cm
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% EI1		float		lower integral limit in energy
% EI2		float		upper integral limit in energy
% filename	string		name of the normal TRANSP output (i.e. not fast ions)
% tr		float		time at which the boundary is requested
% RR		real		Requested radial coordinate to plot the FID
% RZ		real		Requested vertical position to plot the FID	
% plotyn	integer		0 = no plot, 1 = polt
% CM 		real		Z max of the FID
% IE		float		injection energy in eV
%
% OUTPUT
%
% FIDDST		array		trapped fast ion density in the energy range
% FIDDSP		array		passing fast ion density integrated over the specified region
% FIDDSTT		float		trapped total number of fast ions
% FIDDSPT		float		passing total number of fast ions 
% id		integer		array index of FID which is closest to the requested point in the poloidal plane


% Create the array for plotting the 2D FID
[pa, en] = meshgrid(PA, EN);
[X Y] = meshgrid (linspace(0,200,200), linspace(-100, 100, 200));


% Check the integration limits
if (EI1 < min(EN))
  EI1 = min(EN);
endif

% Find the indexes of the integral limits in the array
i1 = max(find(EN <= EI1));
  ENS1 = EN(i1);
i2 = max(find(EN <= EI2));
  ENS2 = EN(i2);

% Gets the trapped - passing boundary and the norm poloidal flux coordinate
plotyn = 0;
[X, BMIN, BMAX, TPB, tb, xb, TPBtb] = fastions_tp_boundary(filename, tr, plotyn);
clear X BMIN BMAX TPB  
  
% Find the required normalized flux surface by determining the minor radius a for a given point
% in the poloidal plane, where a is the distance from the magnetic asis to the LCFS. Once
% the xr = d/a is calculated, the position in the Trapped Passing Boundary is searched and
% the integration ranges in the pitch angles are set symmetrically to zero.
RM = mean(R2D(1:4)); 					% radial coordinate of the magnetic axis
ZM = mean(Z2D(1:4)); 					% vertical coordinate of the magnetic axis
thetab = atan2(Z2D(761:end)-ZM, R2D(761:end)-RM);	% poloidal angle of the last flux surface (boundary)
thetar = atan2(Z2D-ZM, R2D-RM);			% poloidal angles of all the points (requested)
d = sqrt((R2D-RM).^2 + (Z2D-ZM).^2);			
a = sqrt((R2D(761:end)-RM).^2 + (Z2D(761:end)-ZM).^2);
for k = 1:length(thetar)
  it = max(find(thetab <= thetar(k)));
  xr(k) = d(k)/a(it);
  ir = max(find(xb <= xr(k))); 
  pa2(k) = TPBtb(ir);
  pa1(k) = -pa2(k);
  i3(k) = max(find(PA <= pa1(k)));			% pitch angle integration lower limit
  i4(k) = max(find(PA <= pa2(k)));			% pitch angel integration upper limit
end
  
  
% Integration step size
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);

% Calculates the trapped fast ion density in the specified energy/pitch angle range
for k = 1:length(FID)
  w(k) = sum(sum(FID(k,i3(k):i4(k),i1:i2),3),2)*DPA*DEN;
end
w = w';
FIDDST = griddata(R2D, Z2D, w, X, Y, 'linear');
u = find(isnan(FIDDST) == 1);
FIDDST(u) = 0;
clear u w

% Calculates the passing fast ion density in the complement of the specified energy/pitch angle range
for k = 1:length(FID)
  w(k) = sum(sum(FID(k,[1:i3(k)-1, i4(k)+1:end],i1:i2),3),2)*DPA*DEN;
end
w = w';
FIDDSP = griddata(R2D, Z2D, w, X, Y, 'linear');
u = find(isnan(FIDDSP) == 1);
FIDDSP(u) = 0;
clear u w

% Calculates the maximum for the plot
cmax = max(max(max(FIDDST)), max(max(FIDDSP)));


% Calculates the fast ion distribution function limited to the specified energ/pitch angle range
% this is used only for ploting the trapped and passing FID for a given point on the polodal plane 
FIDST = zeros(size(FID));	% trapped
FIDSP = zeros(size(FID));	% passing
for k = 1:length(FID)
  FIDST(k,i3(k):i4(k), i1:i2) = FID(k,i3(k):i4(k), i1:i2);
  FIDSP(k, [1:i3(k)-1, i4(k)+1:end], i1:i2) = FID(k, [1:i3(k)-1, i4(k)+1:end], i1:i2);
end

% Calculates the total number of fast ions
dr = (X(1,2)-X(1,1));
dz = (Y(2,1)-Y(1,1));
FIDDSTT = 2*pi*sum(sum(FIDDST.*X))*dr*dz;
FIDDSPT = 2*pi*sum(sum(FIDDSP.*X))*dr*dz;



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
  figure(3, 'position', [100 100 800 800])
    subplot(2,2,1)
      plot(R2D, Z2D, '.', R2D(id), Z2D(id), 'ro', RR, RZ, 'ks')
      xlabel('R (cm)')
      ylabel('Z (cm)')
      legend('FID points', 'Closest FID', 'Requested FID')
      axis([0 200 -100 100])
      axis equal
      title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])

    subplot(2,2,2)
      pcolor(pa, en, squeeze(FID(id,:,:))');
      shading interp
      hold on
	plot([pa1(id) pa1(id)], [0 IE], '-w', 'linewidth', 2)
	plot([pa2(id) pa2(id)], [0 IE], '-w', 'linewidth', 2)
      hold off
      xlabel('pitch angle')
      ylabel('energy eV')
      axis([-1 1 0 IE 0 CM])
      caxis([0 CM])
      %axis tight
      title(['FID Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.'])
      colorbar      
      
    subplot(2,2,3)
      pcolor(pa, en, squeeze(FIDST(id,:,:))');
      shading interp
      xlabel('pitch angle')
      ylabel('energy eV')
      axis([-1 1 0 IE 0 CM])
      caxis([0 CM])
      %axis tight
      title(['Trapped Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.'])
      colorbar

    subplot(2,2,4)
      pcolor(pa, en, squeeze(FIDSP(id,:,:))');
      shading interp
      xlabel('pitch angle')
      ylabel('energy eV')
      axis([-1 1 0 IE 0 CM])
      caxis([0 CM])
      %axis tight
      title(['Passing Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.'])
      colorbar
      
    
    
load mycolormaps     
  figure(4, 'position', [100 100 800 800])
    subplot(2,2,1)
      pcolor(X, Y, FIDDST);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 cmax])
      title(['Trapped FID - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar
      
    subplot(2,2,2)
      pcolor(X, Y, FIDDSP);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 cmax])
      title(['Passing FID - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar
       
    subplot(2,2,3)
      surf(X, Y, FIDDST);
      shading interp
      %set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      zlabel('Fast Ion Density')
      axis([40 160 -70 70 0 cmax])
      view(10,60)
      caxis([0 cmax])
      colorbar
      title(['Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
    
    subplot(2,2,4)
      surf(X, Y, FIDDSP);
      shading interp
      %set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      zlabel('Fast Ion Density')
      axis([40 160 -70 70 0 cmax])
      view(10,60)
      caxis([0 cmax])
      colorbar
      title(['Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])        
    
endif














	
