function [yntt_na, BBN4_NA, BTN4_NA, yntti_na, YNTA, r_zero, yntt_na_zero] =nd_special_emissivity_JET(r, z,  ynth, R, Z, BBNTX, BTNTX, plotyn, saveyn, filename);
% Function that sum the special non flux averaged neutron emisiivity from BB and BT
% components to the flux averaged thermal neutron emissivity
%
% INPUT
% r		array		time averaged radial coordinate in cm
% z		array		time averaged vertical coordinate in cm
% ynth		array 		time averaged thermal neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% RNA		array		radial coordinates of the non averaged flux surfaces in a given time in cm
% ZNA		array		vertical coordinates of the non averaged flux surfaces in a given time in cm
% BBN4		array		non-flux averaged beam-beam neutron emissivity in cm-3 s-1
% BTN4		array		non-flux averaged beam-thermal neutron emissivity in cm-3 s-1
% plotyn	integer		keyword: 0 no plot, 1 plot
% saveyn 	integer		keyword; 0 no save, 1 save
% filename	string		name of the NETCDF file containing TRANSP output
%
% OUTPUT
% yntt_na	array 		non flux averaged, time averaged total neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% BBN4_NA 	array 		non-flux averaged beam-beam neutron emissivity in cm-3 s-1 interpolated on the regular grid
% BTN4_NA	array		non-flux averaged beam-thermal neutron emissivity in cm-3 s-1 interpolated on the regular grid
% YNTA		float		total neutron yield in cm-3 s-1


% Interpolate the data on the grid used by the standard TRANSP output
BBN4_NA = griddata(R, Z , BBNTX, r, z);
w = find(isnan(BBN4_NA)==1);
BBN4_NA(w) = 0;
BTN4_NA = griddata(R, Z , BTNTX, r, z);
w = find(isnan(BTN4_NA)==1);
BTN4_NA(w) = 0;
% Calculates the total neutron emissivity
yntt_na = BBN4_NA + BTN4_NA + ynth;
u = isnan(yntt_na);
yntt_na(u) = 0;
clear u

% Interpolates on a regular fine mesh and returns the profiles
[ri zi] = meshgrid(20:2:150, -140:2:140);
yntti_na = griddata(r, z, yntt_na, ri, zi);
u = isnan(yntti_na);		% find the indexes in yntti of NAN ...
yntti_na(u) = 0;		% and set them to ZERO
clear u

% Calculates the total neutron yield in the time interval
dr = ri(1,2) - ri(1,1);
dz = zi(2,1) - zi(1,1);
YNTA = 2*pi*sum(sum(yntti_na.*ri))*dr*dz;


% Find the Z = 0 radial profile
r_zero = linspace(0, 200, 400);
yntt_na_zero = interp2(ri, zi, yntti_na, r_zero, 0);
u = isnan(yntt_na_zero);		
yntt_na_zero(u) = 0;			
clear u

return
if (plotyn == 1)
  
  
  load mycolormaps
  figure(3, 'position', [800 200 1000 800])
		pcolor(r,z, yntt_na)
		shading interp
		hold on
		  for k = 1:21;
		     plot(r(k,:), z(k,:), 'k')
		  end
		hold off
		set(gcf,'Colormap',MASTcmap) 
		xlabel('R (cm)','fontweight','bold','fontsize',20)
		ylabel('Z (cm)','fontweight','bold','fontsize',20)
		axis([0 240 -140 140])
		axis equal
		caxis([0 max(max(yntt_na))])
		title('Tota non-flux averaged emissivity cm^{-3} s^{-1}','fontweight','bold','fontsize',20)
		c = colorbar;
    caxis([0 1.5e8])
    set(c, 'fontweight','bold','fontsize',16)
   
    set(gca,'linewidth',4, 'fontsize', 20);
		%text (100, 180, strrep(filename, '_', ' '))
    
    endif
    
   return 
   
% Plot the results. 
if (plotyn == 1)
  load mycolormaps
  figure(30, 'position', [800 200 1000 800])
	subplot(2,2,1)
		pcolor(r,z, yntt_na)
		shading interp
		hold on
		  for k = 1:21;
		     plot(r(k,:), z(k,:), 'k')
		  end
		hold off
		set(gcf,'Colormap',MASTcmap) 
		xlabel('R (cm)')
		ylabel('Z (cm)')
		axis([0 240 -140 140])
		axis equal
		caxis([0 max(max(yntt_na))])
		title('Total non-flux averaged emissivity cm^{-3} s^{-1}')	
		colorbar
		text (100, 180, strrep(filename, '_', ' '))

	subplot(2,2,2)
		pcolor(r,z, BTN4_NA)
		shading interp
		hold on
		  for k = 1:21;
		     plot(r(k,:), z(k,:), 'k')
		  end
		hold off
		set(gcf,'Colormap',MASTcmap) 
		xlabel('R (cm)')
		ylabel('Z (cm)')
		axis([0 240 -140 140])
		axis equal
		caxis([0 max(max(yntt_na))])
		title('Beam Target non-flux averaged emissivity cm-3 s-1')	
		%colorbar

	subplot(2,2,3)
		pcolor(r,z, BBN4_NA)
		shading interp
		hold on
		  for k = 1:21;
		     plot(r(k,:), z(k,:), 'k')
		  end
		hold off
		set(gcf,'Colormap',MASTcmap) 
		xlabel('R (cm)')
		ylabel('Z (cm)')
		axis([0 240 -140 140])
		axis equal
		caxis([0 max(max(yntt_na))])
		title('Beam Beam non-flux averaged emissivity cm-3 s-1')	
		%colorbar

	subplot(2,2,4)
		pcolor(r,z, ynth)
		shading interp
		hold on
		  for k = 1:21;
		     plot(r(k,:), z(k,:), 'k')
		  end
		hold off
		set(gcf,'Colormap',MASTcmap) 
		xlabel('R (cm)')
		ylabel('Z (cm)')
		axis([0 240 -140 140])
		axis equal
		caxis([0 max(max(yntt_na))])
		title('Thermal emissivity flux averaged cm-3 s-1')	
		%colorbar
		
		
  figure(51)
    plot(r_zero, yntt_na_zero, 'r','linewidth', 2)
    xlabel('minor radius (m)')
    ylabel('Total neutron emissivity (cm-3 s-1)')
    legend('Non flux averaged profile')
    title(strrep(filename, '_', ' '));
    
endif

% Save the total neutron emissivity
if (saveyn == 1)
  save('-binary', 'transp_yntt_na.mat', 'r', 'z', 'yntt_na');
  save('-binary', 'transp_yntti_na.mat', 'ri', 'zi', 'yntti_na');
  
  % Radial neutron emissivity
  w = [r_zero' yntt_na_zero'];
  save('-ascii', 'transp_yn_na_radial.dat', 'w');  
endif


