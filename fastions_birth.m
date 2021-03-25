#function [x, y, r, z, theta, xx, yy, FIB, M, R, Z] = fastions_birth(filename, plotyn, saveyn)
% Scripts that reads and makes a 2D histogram of the
% NBI fast ion birth location from TRANSP/NUBEAM runs
close all;
clear all;
% Filename
#filename = '/home/andrea/TRANSP/TRANSP_analysis_script/JET/94665/94665P04_birth.cdf1';
filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_birth.cdf1';
#filename = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29210/O22/29210O23_birth.cdf6';

load('RBW.mat')

pitch = nc_read(filename, 'bs_xksid_D_MCBEAM');
E = nc_read(filename, 'bs_einj_D_MCBEAM');
weightD = nc_read(filename, 'bs_wght_D_MCBEAM');
#weightT = nc_read(filename, 'bs_wght_T_MCBEAM');
%n=40;
% Radial position in (cm) 
r = nc_read(filename, 'bs_r_D_MCBEAM');
rgc = nc_read(filename, 'bs_rgc_D_MCBEAM');
 
% Vertical position in (cm)
z = nc_read(filename, 'bs_z_D_MCBEAM');
%z = repmat(z,n,1);
zgc = nc_read(filename, 'bs_zgc_D_MCBEAM');

figure(1)
scatter(r,z,'r')
hold all
scatter(rgc, zgc, 'b')
xlabel ('R (cm)', 'fontsize', 12);
ylabel ('Z (cm)', 'fontsize', 12);
legend('prt position','gc')

figure(2)
scatter(rgc-r,zgc-z,'r')
legend('gc-prt position')

xlabel ('R (cm)', 'fontsize', 12);
ylabel ('Z (cm)', 'fontsize', 12);

return
% Toroidal angle in deg
theta = nc_read(filename, 'bs_zeta_D_MCBEAM');
%theta = repmat(theta,n,1);

% Calculates the carthesian coordinates
x = r.*cos(theta*2*pi/360);
y = r.*sin(theta*2*pi/360);

M = [x y];

% Generate the histogram
#addpath('/home/andrea/Documents/MAST/Octave/')
edges = linspace(-350, 350, 100);
FIB = hist2d(M,edges, edges);
[R, Z] = meshgrid(edges, edges);

return
figure(3)
    pcolor(R/100, Z/100, FIB); 
    colormap(RBW);
    shading interp;    
    xlabel ('x', 'fontsize', 12);
    ylabel ('y', 'fontsize', 12);
    h = colorbar;
    set(h, 'title', 'Counts')
    title( strrep(filename, '_', ' '), 'fontsize', 12)
    xlim([-4,4])
    ylim([-4,4])
    axis equal

