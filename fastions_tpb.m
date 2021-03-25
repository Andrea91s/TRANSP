function [FIDN, FIDNT, FIDtrapN, FIDtrapNT, FIDtrapR, FIDtrapTR, FIDpassN, FIDpassNT, FIDpassR, FIDpassTR] = fastions_tpb(R2D, Z2D, PA, EN, FID, BMVOL, EI1, EI2, filename, tr, RR, RZ, plotyn, CM, IE)
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
% BMVOL		array		zone volumes
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
% FIDN		array		number of fast ion in each zone
% FIDNT		float		total number of fast ions
% FIDtrapN	array		number of trapped fast ions in each zone
% FIDtrapNT	float		total number of trapped fast ions
% FIDtrapR	array		ratio of trapped to fast ions in each zone
% FIDtrapTR	float		ratio of total trapped to total fast ions
% FIDpassN	array		number of passing fast ions in each zone
% FIDpassNT	float		total number of passing fast ions
% FIDpassR	array		ratio of passing to fast ions in each zone
% FIDpassTR	float		ratio of total passing to total fast ions


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
plotsn = 0;
[XX, BMIN, BMAX, TPB, tb, xb, TPBtb] = fastions_tp_boundary(filename, tr, plotsn);
clear XX BMIN BMAX TPB  
  
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
  i4(k) = max(find(PA <= pa2(k)));			% pitch angle integration upper limit
end

% Integration step size
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);


% Calculates the fast ion distribution function limited to the specified energ/pitch angle range
% this is used only for plotting the trapped and passing FID for a given point on the polodal plane 
FIDST = zeros(size(FID));	% trapped
FIDSP = zeros(size(FID));	% passing
for k = 1:length(FID)
  FIDST(k,i3(k):i4(k), i1:i2) = FID(k,i3(k):i4(k), i1:i2);
  FIDSP(k, [1:i3(k)-1, i4(k)+1:end], i1:i2) = FID(k, [1:i3(k)-1, i4(k)+1:end], i1:i2);
end


% Calculate the total number of fast ions in each zone
% and the total number summing over all the zones
FIDN = zeros(length(FID),1);
%FIDN = sum((sum(FID,2)*DPA)*DEN,3).*BMVOL;
for k = 1:length(FID)
  FIDN(k) = sum((sum(FID(k,1:end,i1:i2),2)*DPA)*DEN,3).*BMVOL(k);
end
%FID = FID';
FIDNT = sum(FIDN);		

% Calculates the numner of TRAPPED fast ions in each zone
% and the total number summing over all the zones
FIDtrapN = zeros(length(FID),1);
for k = 1:length(FID)
  FIDtrapN(k) = sum((sum(FID(k,i3(k):i4(k),i1:i2),2)*DPA)*DEN,3)*BMVOL(k);
end
%FIDtrapN = FIDtrapN';
FIDtrapNT = sum(FIDtrapN);

% Calculates the total numner of PASSIG fast ions in each zone
% and the total number summing over all the zones
FIDpassN = zeros(length(FID),1);
for k = 1:length(FID)
  FIDpassN(k) = sum((sum(FID(k,[1:i3(k)-1, i4(k)+1:end],i1:i2),2)*DPA)*DEN,3).*BMVOL(k);
end
%FIDpassN = FIDpassN';
FIDpassNT = sum(FIDpassN);


% Calculates the ratio of trapped and passing particles in each zone
FIDtrapR = FIDtrapN./FIDN;
FIDpassR = FIDpassN./FIDN;
FIDtrapTR = FIDtrapNT/FIDNT;
FIDpassTR = FIDpassNT/FIDNT;


% Output the results
printf('Total number of fast ions in the range %f to %f keV:\n', EI1/1000, EI2/1000);
printf('  Total: %1.4g\n', FIDNT);
printf('  Trapped: %1.4g \t trapped fraction = %f\n', FIDtrapNT, FIDtrapNT/FIDNT);
printf('  Passing: %1.4g \t passing fraction = %f\n', FIDpassNT, FIDpassNT/FIDNT);
printf('  Trapped + Passing: %1.4g\n', FIDtrapNT + FIDpassNT);
printf('  Trapped: %1.4g \t trapped fraction = %f\n', FIDtrapNT, FIDtrapTR);
printf('  Passing: %1.4g \t passing fraction = %f\n', FIDpassNT, FIDpassTR);
printf('  Trapped + Passing: %1.4g and in fraction %f\n', FIDtrapNT + FIDpassNT, FIDtrapTR + FIDpassTR);


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

% Print the total, trapped and passing poulation at the selected point
% on the poloidal plane
printf('Number of fast ions at R = %f and Z = %f for E = %f to %f keV:\n', R2D(id), Z2D(id), EI1/1000, EI2/1000);
printf('  Total: %1.4g\n', FIDtrapN(id) + FIDpassN(id));
printf('  Trapped: %1.4g \t trapped fraction = %f\n',  FIDtrapN(id), FIDtrapN(id)/(FIDtrapN(id) + FIDpassN(id)));
printf('  Passing: %1.4g \t passing fraction = %f\n', FIDpassN(id), FIDpassN(id)/(FIDtrapN(id) + FIDpassN(id)));
printf('  Total: %1.4g and in fraction %f\n', FIDN(id), FIDtrapN(id)/(FIDtrapN(id) + FIDpassN(id)) + FIDpassN(id)/(FIDtrapN(id) + FIDpassN(id)));
printf('  Trapped: %1.4g \t trapped fraction = %f\n',  FIDtrapN(id), FIDtrapR(id));
printf('  Passing: %1.4g \t passing fraction = %f\n', FIDpassN(id), FIDpassR(id));
printf('  Trapped + Passing: %1.4g and in fraction %f\n', FIDtrapN(id) + FIDpassN(id), FIDtrapR(id) + FIDpassR(id));



% Plot the time traces
if (plotyn == 1)

  [pa, en] = meshgrid(PA, EN);

  % Plot the TRANSP grid and the fast ion distribution in the
  % requested point on the poloidal cross-section
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
      
    
    
  % Plot the trapped and passing fast ion distribution on
  % the poloidal cross section and the radial profile
  % along the major radius
  
  % Create the array for plotting the 2D FID
  [X Y] = meshgrid (linspace(0,200,200), linspace(-100, 100, 200));

  % Interpolates the fast iond distribution on the grid
  FIDtrap = griddata(R2D, Z2D, FIDtrapN./BMVOL, X, Y, 'linear');
      u = find(isnan(FIDtrap) == 1);
      FIDtrap(u) = 0;
      clear u
      
  FIDpass = griddata(R2D, Z2D, FIDpassN./BMVOL, X, Y, 'linear');
      u = find(isnan(FIDpass) == 1);
      FIDpass(u) = 0;
      clear u      
  
  FIDtrap_ratio = griddata(R2D, Z2D, FIDtrapR, X, Y, 'linear');
      u = find(isnan(FIDtrap_ratio) == 1);
      FIDtrap_ratio(u) = 0;
      clear u
      
  FIDpass_ratio = griddata(R2D, Z2D, FIDpassR, X, Y, 'linear');
      u = find(isnan(FIDpass_ratio) == 1);
      FIDpass_ratio(u) = 0;
      clear u

  FIDtotal = griddata(R2D, Z2D, FIDN./BMVOL, X, Y, 'linear');
      u = find(isnan(FIDtotal) == 1);
      FIDtotal(u) = 0;
      clear u        
      
  % Calculates the maximum for the plot
  cmax = max(max(max(FIDtrap)), max(max(FIDpass)));

  
  figure(4, 'position', [100 100 800 800])
  load mycolormaps   
    subplot(2,2,1)
      pcolor(X, Y, FIDtrap);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 cmax])
      title(['Trapped FI - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar
      
    subplot(2,2,2)
      pcolor(X, Y, FIDpass);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 cmax])
      title(['Passing FI - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar
       
    subplot(2,2,3)
      pcolor(X, Y, FIDtrap_ratio);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 1])
      title(['Ratio of Trapped FI - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar
      
    subplot(2,2,4)
      pcolor(X, Y, FIDpass_ratio);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -70 70])
      axis equal
      caxis([0 1])
      title(['Ratio of Passing FI - Energy Range ' num2str(ENS1/1000) ' to ' num2str(ENS2/1000) ' keV.' ])
      colorbar   

      
  % Plot of the radial profiles    
      
  % Extract the radial profile of the trapped, passing and total fast ions
  r = linspace(0, 200, 400);
  z = zeros(1,400);
  fid_trap = interp2(X, Y, FIDtrap, r, z);
  fid_pass = interp2(X, Y, FIDpass, r, z);
  fid_total = interp2(X, Y, FIDtotal, r, z);
  trap_ratio = fid_trap./fid_total;
    u = find(isinf(trap_ratio) == 1);
    trap_ratio(u) = 0;
  pass_ratio = fid_pass./fid_total;
    u = find(isinf(pass_ratio) == 1);
    pass_ratio(u) = 0;
  
  figure(5, 'position', [100 100 800 400])
    subplot(1,2,1)
      plot(r, fid_trap, 'r', 'linewidth', 2, r, fid_pass, 'b', 'linewidth', 2, r, fid_pass + fid_trap, 'm', 'linewidth', 2, r, fid_total, 'k', 'linewidth', 2)
      xlabel('major radius (cm)')
      ylabel('Fast Ions')
      legend('TRAPPED', 'PASSING', 'PASS + TRAP', 'TOTAL')
    subplot(1,2,2)
      plot(r, trap_ratio, 'r', 'linewidth', 2, r, pass_ratio, 'b', 'linewidth',2)
      xlabel('major radius (cm)')
      ylabel('Fast Ions Fractions')
      legend('TRAPPED', 'PASSING')
endif














	
