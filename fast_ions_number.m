function [FI] = fast_ions_number(filename, ts, plotyn, saveyn)
% Function that calculates the beam ions density from the standard
% NetCDF file using the BDENS (1/cm3) and DVOL (cm3) quantities
%
% Example:
% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18.CDF';
% fast_ions_number(filename, ts = 0.216, plotyn = 1, saveyn = 0);

MNEUT = nc_read(filename, 'MNEUT');  %MEASURED NEUTRONS
BDENS = nc_read(filename, 'BDENS');  % BEAM ION DENSITY
BMVOL = nc_read(filename, 'BMVOL');  % VOLUMES volume
DVOL = nc_read(filename, 'DVOL');   % ZONE VOLUME
TIME = nc_read(filename, 'TIME');
PVOL = nc_read(filename, 'PVOL');  % PLASMA VOLUME
X = nc_read(filename, 'X');
TRAPB_D = nc_read(filename, 'TRAPB_D');

% Calculates the FI number
FIN = sum(BDENS.*DVOL,2);
FIND = BDENS.*DVOL;

% Calcualtes the FI density
FID = FIN./PVOL;

% Extract the FI density profile for the selected time
idx = max(find(TIME <= ts));
idx=56;

x = X(idx,:);
fid = BDENS(idx,:);

% Copy everything into the return structure
FI.filename = filename;
FI.MNEUT = MNEUT;
FI.selected_time = ts;
FI.transp_time = TIME(idx);
FI.x = x;
FI.fid = fid;
FI.TIME = TIME;
FI.BDENS = BDENS;
FI.FIN = FIN;
FI.FID = FID;
FI.TRAPB_D = TRAPB_D;

% plot
if (plotyn == 1)
figure(20)
plot(TIME, MNEUT, 'linewidth', 2)
    xlabel ('time (s)')
    ylabel ('Neutrons yield')
   % legend('\Sigma_j BDENS_{i,j} \times DVOL_{i,j}')
    title(filename)
    
figure(10, 'position', [100 100 560 630])
  subplot(3,1,1)
    plot(TIME, FIN, 'linewidth', 2)
    xlabel ('time (s)')
    ylabel ('number of fast ions')
    legend('\Sigma_j BDENS_{i,j} \times DVOL_{i,j}')
    title(filename)
  subplot(3,1,2)
    plot(TIME, FID, 'linewidth', 2)
    xlabel ('time (s)')
    ylabel ('fast ions density (cm^{-3})')
    legend('(\Sigma_j BDENS_{i,j} \times DVOL_{i,j})/PVOL_i')
    title(filename)    
  subplot(3,1,3)
    plot(TIME, TRAPB_D, 'linewidth', 2)
    xlabel ('time (s)')
    ylabel ('Trapped Fraction')
    legend('TRAPB_D', 'location', 'southeast')
    title(filename) 

   
    
figure(11)
  [m, n] = size(X);
  TT = repmat(TIME, 1, n); 
  pcolor(TT,X,BDENS)
  shading interp
  xlabel('time (s)')
  ylabel('normalized pol. flux')
  title(filename)  
  h = colorbar;
  set(h, 'title', 'BDENS D (cm^{-3})')
  axis tight
 

figure(12)
  pcolor(TT,X,FIND)
  shading interp
  xlabel('time (s)')
  ylabel('normalized pol. flux')
  h = colorbar;
  set(h, 'title', '#')
  axis tight
  
figure(13)
  plot(x, fid, 'k', 'linewidth', 2);
  xlabel('normalized pol. flux')
  ylabel ('fast ions density (cm^{-3})')
  title(filename)        
  legend(['FI profile for t = ' num2str(TIME(idx)) ' s.'])
  
endif