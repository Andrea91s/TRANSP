function [yn_asym] = nd_emissivity_asymmetry(r, z,  yntt, yntt_na, ri, zi, yntti, yntti_na, plotyn, saveyn)
% Function that plots the difference between the flux averaged and non flux averaged neutron emissivity

% Calculates the difference
yn_asym = yntt_na - yntt;
yni_asym = yntti_na - yntti;

% Plot the results. 
if (plotyn == 1)
  load mycolormaps
  figure(4)
      yy1 = max([max(max(yn_asym)), abs(min(min(yn_asym)))]);
      load DIFFmap.mat
      pcolor(r,z, yn_asym)
      shading interp
      hold on
	for k = 1:21;
	    plot(r(k,:), z(k,:), 'k')
	end
      hold off
      set(gcf,'Colormap',DIFFmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([0 240 -140 140])
      axis equal
      caxis([-yy1 yy1])
      title('Difference averaged emissivity cm-3 s-1')	
      colorbar

  figure(40)
      yy1 = max([max(max(yni_asym)), abs(min(min(yni_asym)))]);
      load DIFFmap.mat
      pcolor(ri,zi, yni_asym)
      shading interp
      hold on
	for k = 1:21;
	    plot(r(k,:), z(k,:), 'k')
	end
      hold off
      set(gcf,'Colormap',DIFFmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([0 240 -140 140])
      axis equal
      caxis([-yy1 yy1])
      title('Difference averaged emissivity cm-3 s-1')	
      colorbar      
endif

% Save the total neutron emissivity
if (saveyn == 1)
  save('-binary', 'transp_yn_asym.mat', 'r', 'z', 'yn_asym');
endif




