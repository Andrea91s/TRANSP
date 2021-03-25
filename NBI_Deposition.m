filename_global = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18.CDF';
filename_fi = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_fi_1.cdf';
filename_fib = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_birth.cdf1';
% Fake xtor 
x = linspace(0,1,60);

% Time [s]
% 376x1
TIME = nc_read(filename_global, 'TIME');

% Norm. coordinate
% 376x60
X = nc_read(filename_global , 'X');

% FAST ION SCE:  DEPOSITION [1/s]
% 376x1
SFDEP = nc_read (filename_global , 'SFDEP');

% Volumes [cm3]
% 376x60
DV = nc_read (filename_global , 'DVOL');

% D BEAM DEPOSITION (TOTAL) [1/cm3 1/s]
% 376x60
BDEP_D = nc_read (filename_global , 'BDEP_D');

% Total number of D 
NTOTFI = nc_read (filename_fi, 'NTOT_D_NBI');
BDENS2 = nc_read (filename_fi, 'bdens2');

% Beam ion density [1/cm3]
% 376x60
BDENS = nc_read (filename_global, 'BDENS');

% Beam deposition in [1/cm3 1/s]
% 376x60
BD1 = nc_read(filename_global , 'BDEP01_TOT');
%BD2 = nc_read(filename_global , 'BDEP02_TOT');
%BD3 = nc_read(filename_global , 'BDEP03_TOT');
%BD4 = nc_read(filename_global , 'BDEP04_TOT');
BD = BD1;% + BD2 + BD3 + BD4;

% D BEAM TOTAL DEPOSITION SCE [1/s]
% 376x1
SBDEPSC_D = nc_read (filename_global , 'SBDEPSC_D');

% Beam deposition rate [1/s]
% 376x1
BDR = sum(BD.*DV,2);

% Beam deposition density [1/m3]
% 1x60
BDD = mean(BD.*DV,1);

% D BEAM DEPOSITION: BEAM-BEAM Impact Ionization [1/cm3 1/s]
% 376x60
SDBBI_D = nc_read (filename_global, 'SDBBI_D');

% D BEAM DEPOSITION: BEAM-BEAM CX [1/cm3 1/s]
% 376x60
SDBBX_D = nc_read (filename_global, 'SDBBX_D');

% D BEAM DEP: IONIZ. on electrons [1/cm3 1/s]
% 376x60
SDBIE_D = nc_read (filename_global, 'SDBIE_D');

% D BEAM DEP: IONIZ. on therm.ions [1/cm3 1/s]
% 376x60
SDBII_D = nc_read (filename_global, 'SDBII_D');

% D BEAM DEP: IONIZ. on impurities [1/cm3 1/s]
% 376x60
SDBIZ_D = nc_read (filename_global, 'SDBIZ_D');

% D BEAM DEPOSITION: CX W/D + IONS [1/cm3 1/s]
% 376x60
SDCXD_D = nc_read (filename_global, 'SDCXD_D');

% Adding all the DEPOSITION terms above [1/cm3 1/s]
% This is equal to BD
% 376x60
SDEPTOT = SDBBI_D + SDBBX_D + SDBIE_D + SDBII_D + SDBIZ_D + SDCXD_D;

% Reading the slowing down time
% 376x60
TSL1A_D = nc_read (filename_global , 'TSL1A_D');

% Beam ion density calculated from the Deposition [1/cm3]
% 376x60
FID = SDEPTOT.*TSL1A_D;

% Read the weights from the birth file  [-]
w = nc_read(filename_fib, 'bs_wght_D_MCBEAM');

% Time interval used to dump the FI distribution [1/s]
DT = 1E-3;

% Estimates the deposition rate
EDR = (sum(w)/DT)
NTOT_D_NBI = ncread(filename_fi, 'NTOT_D_NBI')
TSL1A_D(56,:).*DVOL(56,:)
% ----------------------------------------------------------------------------------------------------

% Figures

% Deposition rate
figure
plot (TIME, SFDEP, TIME, BDR, TIME, sum(BDEP_D.*DV,2), TIME, sum(SDEPTOT.*DV,2))
xlabel('time (s)', 'fontsize', 12)
ylabel('rate (1/s)', 'fontsize', 12)
title('Source Deposition Rate')
legend('SFDEP', 'BDR', 'from BDEP\_D', 'from SDEPTOT')

% Deposition rate density profile components
figure 
semilogy(x, mean(SDBBI_D,1), ...
     x, mean(SDBBX_D,1), ...
     x, mean(SDBIE_D,1), ...
     x, mean(SDBII_D,1), ...
     x, mean(SDBIZ_D,1), ...
     x, mean(SDCXD_D,1), ...
     x, mean(SDEPTOT,1), ...
     x, mean(BDEP_D,1), 'o')
xlabel('norm. tor. coord ', 'fontsize', 12)
ylabel('density rate (1/s 1/cm3)', 'fontsize', 12)
title('Source Deposition Density Rate Components')
legend('SDBBI_D', 'SDBBX_D', 'SDBIE_D', 'SDBII_D', 'SDBIZ_D', 'SDCXD_D', 'SDEPTOT', 'BDEP_D')


% Comparing the ion density and the beam deposition
figure
plot (x, mean(BDENS,1), x, mean(FID,1))
xlabel('norm. tor. coord ', 'fontsize', 12)
ylabel('density (1/cm3)', 'fontsize', 12)
title('Beam Ion Density ')
legend('BDENS', 'BDEP\_D x TSL1A\_D')































