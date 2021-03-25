filename = '/home/andrea/TRANSP/AGC_fi_1.cdf';

RR = 70;
RZ = -20;
CM = 0;
IE = 80000;
plotyn=1;
saveyn=0;

close all;

% Read the pitch angle scattering values
PA = nc_read(filename, 'A_D_NBI');

% Read the energies
EN = nc_read(filename, 'E_D_NBI');

% Read the fast ion distribution function
FID = nc_read(filename, 'F_D_NBI');

FID_GLOBAL = squeeze(sum(FID,1));

% Read the R and Z coordinates
R2D = nc_read(filename, 'R2D');
Z2D = nc_read(filename, 'Z2D');


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



% Calculates the fast ion distribution versus energy integrating over
% all pitch angles
FIDE = squeeze(sum(FID,2))*DPA;
FIDP = squeeze(sum(FID,3))*DEN;


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



  figure(1, 'position', [100 100 800 800])
    subplot(2,2,1)
      plot(R2D, Z2D, '.', R2D(id), Z2D(id), 'ro', RR, RZ, 'ks')
      xlabel('R (cm)')
      ylabel('Z (cm)')
      legend('FID points', 'Closest FID', 'Requested FID')
      axis([0 200 -150 150])
      axis equal
      title(filename)

    subplot(2,2,1)
    %{
      pcolor(pa, en, squeeze(FID(id,:,:))');
      shading interp
      xlabel('pitch')
      ylabel('energy eV')
      axis([-1 1 0 IE 0 CM])
      %caxis([0 CM])
      title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
      colorbar
     %} 
   subplot(2,2,3) 
    plot(EN/1000, FIDE(id,:), 'linewidth', 2,'r')
    xlabel ('Energy (keV)')
    ylabel ('Normalized dFID(E)/dE density')
    axis([0 IE/1000])
    #axis([0 IE/1000 0 1.1*max(FIDE(id,:))])
    title(['FID at R = ' num2str(RR) ' cm, Z = ' num2str(RZ) ' cm.'])     
    
    
  subplot(2,2,4) 
    plot(PA, FIDP(id,:), 'linewidth', 2,'r')
    xlabel ('pitch')
    ylabel ('Normalized dFID(Pa)/dPa density')
    axis([-1 1])
    #axis([-1 1 0 1.1*max(FIDP(id,:))])
    title(['FID at R = ' num2str(RR) ' cm, Z = ' num2str(RZ) ' cm.'])       

    
    
    
filename1 = '/home/andrea/TRANSP/AGO_fi_1.cdf';
% Read the pitch angle scattering values
PA = nc_read(filename1, 'A_D_NBI');
% Read the energies
EN = nc_read(filename1, 'E_D_NBI');
% Read the fast ion distribution function
FID = nc_read(filename1, 'F_D_NBI');
FID_GLOBAL = squeeze(sum(FID,1));
% Read the R and Z coordinates
R2D = nc_read(filename1, 'R2D');
Z2D = nc_read(filename1, 'Z2D');
% Read the volume of each zone
BMVOL = nc_read(filename1, 'BMVOL');
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
% Calculates the fast ion distribution versus energy integrating over
% all pitch angles
FIDE = squeeze(sum(FID,2))*DPA;
FIDP = squeeze(sum(FID,3))*DEN;
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
hold all;    
    
    
figure(1, 'position', [100 100 800 800])
subplot(2,2,2)
      plot(R2D, Z2D, '.', R2D(id), Z2D(id), 'ro', RR, RZ, 'ks')
      xlabel('R (cm)')
      ylabel('Z (cm)')
      legend('FID points', 'Closest FID', 'Requested FID')
      axis([0 200 -150 150])
      axis equal
      title(filename1)
subplot(2,2,3) 
   hold all;
    plot(EN/1000, FIDE(id,:), 'linewidth', 2,'b')
    legend('ASCOT GC','ASCOT GO')
subplot(2,2,4) 
  hold all;
    plot(PA, FIDP(id,:), 'linewidth', 2,'b')
    legend('ASCOT GC','ASCOT GO')

    
    
    


  












	
