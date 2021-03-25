function [birthC27] = fastions_birth(filename, plotyn, saveyn)
% Scripts that reads and makes a 2D histogram of the
% NBI fast ion birth location from TRANSP/NUBEAM runs
 
% Filename
% filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29881/29881O07_birth.cdf1';


birthC27.r = nc_read(filename, 'bs_r_D_MCBEAM');
birthC27.z = nc_read(filename, 'bs_z_D_MCBEAM'); 
birthC27.rgc = nc_read(filename, 'bs_rgc_D_MCBEAM');
birthC27.zgc = nc_read(filename, 'bs_zgc_D_MCBEAM');
birthC27.pitch = nc_read(filename, 'bs_xksid_D_MCBEAM');
birthC27.E = nc_read(filename, 'bs_einj_D_MCBEAM');
birthC27.weight = nc_read(filename, 'bs_wght_D_MCBEAM');
birthC27.phi = nc_read(filename, 'bs_zeta_D_MCBEAM');
birthC27.time  = nc_read(filename, 'bs_time_D_MCBEAM');
birthC27.beamid = nc_read(filename, 'bs_ib_D_MCBEAM');




return


% Vertical position in (cm)
z = nc_read(filename, 'bs_zgc_D_MCBEAM');
#z = repmat(z,n,1);

% Toroidal angle in deg
theta = nc_read(filename, 'bs_zeta_D_MCBEAM');
#theta = repmat(theta,n,1);

% Calculates the carthesian coordinates
x = r.*cos(theta*2*pi/360);
y = r.*sin(theta*2*pi/360);

M= [x y];

% Generate the histogram
addpath('/home/andrea/Documents/MAST/Octave/')
edges = linspace(-200, 200, 100);
FIB = hist2d(M,edges, edges);
[R, Z] = meshgrid(edges, edges);

% plot
if (plotyn == 1)

    pcolor(R./100, Z./100,FIB); 
    colormap(RBW);
    shading interp;    
    xlabel ('R (m)', 'fontsize', 12);
    ylabel ('Z (m)', 'fontsize', 12);
    h = colorbar;
    set(h, 'title', 'Counts')
    title( strrep(filename, '_', ' '), 'fontsize', 12)
    axis equal
endif

% save
if (saveyn == 1)
  r=r./100;
  z=z./100;
save('-ascii', 'R_ascot.txt', 'r');
save('-ascii', 'Z_ascot.txt', 'z');
save('-ascii', 'phi_ascot.txt', 'theta');
endif
