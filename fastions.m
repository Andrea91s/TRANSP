function [R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id] = fastions(filename, RR, RZ, CM, IE, plotyn, saveyn)
% Script to access neutrons data from TRANSP Fast Ion output file
% INPUT
% filename	string		name of the NETCDF file containing TRANSP output
% RR		real		Requested radial coordinate to plot the FID
% RZ		real		Requested vertical position to plot the FID	
% CM 		real		Z max of the FID: 0 to auto
% IE		real		injextion energy in eV
% plotyn 	integer		keyword: 0 no plot, 1 plot
% saveyn	integer		0 does not save the FID, 1 save the FID
%
% OUTPUT
% R2D		array		radial position at which the FID has been calculated in cm
% Z2D		array		vertical position at which the FID has been calculated in cm
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% BMVOL		array		zone volumes
% X		array		grid of radial coordinates for the 2D plots
% Y		array		grid of vertical coordinates for the 2D plots
% FIDD		array		fast ion distribution density
% FIDDT		float		total number of fast ions on the interpolated grid
% FIDN		float		total number of fast ions from FBM
% BDENS2T	float 		total number of fast ions from BDENS2	
% FIDE		array		fast ion distribution function integrated over all pitch angles as a function of energy
% FIDP		array		fast ion distribution function integrated over all energy as a function of pitch ange
% r		array		radial coordinate for Z = 0
% er		array		fid density profile for Z = 0 
% id		integer		array index of FID which is closest to the requested point in the poloidal plane
%
% Example
% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P02/29909P02_fi_1.cdf';
% RR = 100;
% RZ = 0;
% CM = 0;
% IE = 80000;
% plotyn = 1;
% saveyn = 0;
%  [R2D, Z2D, PA, EN, FID, X, Y, FIDD, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id] = fastions(filename, RR, RZ, CM, IE, plotyn, saveyn);


% Read the pitch angle scattering values
PA = nc_read(filename, 'A_D_NBI');
PAB = nc_read(filename, 'AB_D_NBI');

% Read the energies
EN = nc_read(filename, 'E_D_NBI');
ENB = nc_read(filename, 'EB_D_NBI');

% Read the fast ion distribution function
FID = nc_read(filename, 'F_D_NBI');

FID_GLOBAL = squeeze(sum(FID,1));


BDENS2 = nc_read(filename, 'bdens2');
BDENS2 = BDENS2';

% Read the R and Z coordinates
R2D = nc_read(filename, 'R2D');
Z2D = nc_read(filename, 'Z2D');

% Read the volume of each zone
BMVOL = nc_read(filename, 'BMVOL');

% Create the array for plotting the 2D FID
[pa, en] = meshgrid(PA, EN);
[X Y] = meshgrid (linspace(0,200,200), linspace(-100, 100, 200));


% Interpolates the FID on the grid for plotting
DEN = EN(2) - EN(1);
DPA = PA(2) - PA(1);
w = sum(sum(FID,3),2)*DPA*DEN;
FIDD = griddata(R2D, Z2D, w, X, Y, 'linear');
u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u


% Calculates the total number of fast ions in each zone and then
% suns all the zones to get the total number of fast ions
DEN = EN(2) - EN(1);	% energy interval
DPA = PA(2) - PA(1);	% pith angle inteval

FIDN = sum(sum((sum(FID,2)*DPA)*DEN,3).*BMVOL,1)/2;
BDENS2T = sum(BDENS2);
BDENS2T = sum(BDENS2);
dr = (X(1,2)-X(1,1));
dz = (Y(2,1)-Y(1,1));
FIDDT = 2*pi*sum(sum(FIDD.*X))*dr*dz;
printf('Total number of FI from FBM: %1.4g\n', FIDN);
printf('Total number of FI from BDENS2: %1.4g\n', BDENS2T);
printf('Total number of FI from interpolated FID: %1.4g\n', FIDDT);



% Calculates the fast ion distribution versus energy integrating over
% all pitch angles
FIDE = squeeze(sum(FID,2))*DPA;
FIDP = squeeze(sum(FID,3))*DEN;


% Calculate the FID in the vpar, vper space using the boundary
% values and interpolating the FID on them
v = sqrt(ENB);
theta = acos(PAB);
vper = zeros(length(theta),length(v));
vpar = zeros(length(theta),length(v));
for m = 1:length(theta)
  for n = 1:length(v)
    vpar(m,n) = v(n)*cos(theta(m));
  end
end
for m = 1:length(theta)
  for n = 1:length(v)
    vper(m,n) = v(n)*sin(theta(m)); 
  end 
end
vpar = vpar/max(max(vpar));
vper = vper/max(max(vper));
[u, w] = meshgrid(EN, PA);
[U, W] = meshgrid(ENB, PAB);
NF = size(FID)(1);
FIDB = zeros(NF, length(theta), length(v));

for k = 1:NF
  FIDB(k,:,:) = interp2(u,w,squeeze(FID(k,:,:)), U, W);
end

www = find(isnan(FIDB) == 1);
FIDB(www) = 0;
clear www


% Extract the radial profile of the neutron emissivity
r = linspace(0, 200, 400);
z = zeros(1,400);
er = interp2(X, Y, FIDD, r, z);


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
 
    FIDTR = squeeze(FID(id,:,:))';

    #save('ENERGY_TRANSP.dat','EN','-ascii')
    #save('PITCH_TRANSP.dat','PA','-ascii')
    #save('TRANSP_R130_Z0_29981_0245.dat','FIDTR','-ascii')
   # save('TRANSP_GLOBAL_29909_0216_CXOFFP05.dat','FID_GLOBAL','-ascii')
  
% Plot the time traces
if (plotyn == 1)
  figure(1, 'position', [100 100 800 800])
    subplot(2,2,1)
      plot(R2D, Z2D, '.', R2D(id), Z2D(id), 'ro', RR, RZ, 'ks')
      xlabel('R (cm)')
      ylabel('Z (cm)')
      legend('FID points', 'Closest FID', 'Requested FID')
      axis([0 200 -150 150])
      axis equal
      title(filename)

    subplot(2,2,2)
      pcolor(pa, en, squeeze(FID(id,:,:))');
      shading interp
      xlabel('pitch angle')
      ylabel('energy eV')
      axis([-1 1 0 IE 0 CM])
      %caxis([0 CM])
      title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
      colorbar
      
   subplot(2,2,3) 
    plot(EN/1000, FIDE(id,:), 'linewidth', 2,'r')
    xlabel ('Energy (keV)')
    ylabel ('dFID(E)/dE density')
    axis([0 IE/1000])
    #axis([0 IE/1000 0 1.1*max(FIDE(id,:))])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])    
  
  subplot(2,2,4) 
    plot(PA, FIDP(id,:), 'linewidth', 2,'r')
    xlabel ('pitch angle')
    ylabel ('dFID(Pa)/dPa density')
    axis([-1 1])
    #axis([-1 1 0 1.1*max(FIDP(id,:))])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])   
    
   endif 

 
load mycolormaps     
  figure(2, 'position', [100 100 800 400])
    subplot(1,2,1)
      pcolor(X, Y, FIDD);
      shading interp
      set(gcf,'Colormap',MASTcmap) 
      xlabel('R (cm)')
      ylabel('Z (cm)')
      axis([40 160 -100 100])
      axis equal
      title(['Fast Ion Distribution'])
      colorbar
       

    subplot(1,2,2)
      plot(r,er, 'b', 'linewidth',2)
      xlabel('R (cm)')
      ylabel('Fast Ion Distribution Function Density')
  
  figure(7)
      %pcolor(vpar, vper , squeeze(FIDB(id,:,:)));
      contourf(vpar, vper , squeeze(FIDB(id,:,:)), 'linewidth', 2);
      shading interp
      xlabel('Norm. Vpar')
      ylabel('Norm. Vper')
      axis([-1 1 0 1 0 CM])
      axis equal
      caxis([0 CM])
      title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
      colorbar
     
  


  












	
