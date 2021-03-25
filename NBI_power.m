function [TIME, PINJ, BPOH, BPLIM, BPSHI, BPCXI, BPCXE, PBPCX, BPTE, BPTi, BPTH, BPBAL] = NBI_power(filename, plotyn)
% Script to access NBI power balance data from TRANSP output file
% INPUT
% filename	string		name of the NETCDF file containing TRANSP output
% plotyn 	integer		keyword: 0 no plot, 1 plot
%
% OUTPUT
% TIME		array		time in sec
% PINJ      array       injected power [Watts]
% BPOH      array       OH power to fast ions [Watts]
% BPBAL     array       power balance [Watts]
% BPLIM     array       lost orbit power losses [Watts]
% BPSHI     array       shine-through power losses [Watts]
% BPCXI     array       power to CX (inside LCFS) [Watts]
% BPCXX     array       power to CX (outside LCFS) [Watts]
% BPTE      array       power to electrons [Watts]
% BPTI      array       power to ions [Watts]
% BPTH      array       fast ions thermalized power [Watts]
% BPCI0     array       D BEAM CX SCE POWER (INT) WATTS
% BPCOL     array       D BEAM PWR: COLLISIONAL TORQUE WATTS
% BPCPR     array       D POWER: COMPRESSION OF D BEAM WATTS
% BPCRI     array       D BEAM CX RECAPTURE (INT) WATTS
% BPCRX     array       D BEAM CX RECAPTURE (EXT) WATTS
% BPCX0     array       D BEAM CX SCE POWER (EXT) WATTS
% BPST      array       D BEAM POWER STORED WATTS
% 
% Example:
% filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U16/29880U16.CDF';
% NBI_power(filename, plotyn = 1)



% Read the time
TIME = nc_read(filename, 'TIME');   % time in sec

% Input Power: note that some inputs powers are not included (RF,...)
PINJ = nc_read(filename, 'PINJ');
BPOH = nc_read(filename, 'BPOH_D');
BPCPR = nc_read(filename, 'BPCPR_D');

% Power losses 
BPLIM = nc_read(filename, 'BPLIM_D');
BPSHI = nc_read(filename, 'BPSHI_D');
#BPCOL = nc_read(filename, 'BPCOL_D');

% Power losses to CX
BPCXI = nc_read(filename, 'BPCXI_D');
BPCXX = nc_read(filename, 'BPCXX_D');
BPCI0 = nc_read(filename, 'BPCI0_D');
BPCRI = nc_read(filename, 'BPCRI_D');
BPCRX = nc_read(filename, 'BPCRX_D');
BPCX0 = nc_read(filename, 'BPCX0_D');
BPCX = BPCXX + BPCXI + BPCRI + BPCRX + BPCX0 + BPCI0;

% power stored
BPST = nc_read(filename, 'BPST_D');

% Power to thermal electrons and ions
BPTE = nc_read(filename, 'BPTE_D');
BPTI = nc_read(filename, 'BPTI_D');
BPTP = BPTE + BPTI;

% What is left becomes thermalized fast ion powers
BPTH = nc_read(filename, 'BPTH');

% Total power losses
BPLT = BPLIM + BPSHI + BPCX;


% Total power input
BPIT = PINJ + BPOH + BPCPR;

% Power balance check
BPBAL = nc_read(filename, 'BPBAL_D');
#BPbal = BPIT - BPLT - BPCOL - BPTP;



% Plot the time traces
if (plotyn == 1)
  figure(1)
    plot(TIME, BPIT, 'k', 'linewidth', 2, TIME, BPLT, 'r', 'linewidth', 2, TIME, BPTE, 'b', 'linewidth', 2,TIME, BPTI, 'g', 'linewidth', 2)
    axis([TIME(1) TIME(end)])
    xlabel('time (s)')
    ylabel('NBI Power (Watts)')
    legend('Input', 'Total losses', 'To electrons','To ions', 'location', 'northeast') 
    legend('boxoff')
    title(filename)
%{
  figure(2)
    plot(TIME, BPbal, 'k', 'linewidth', 2, TIME, BPBAL, 'r', 'linewidth', 2) 
    axis([TIME(1) TIME(end)])
    xlabel('time (s)')
    ylabel('NBI Power (Watts)')
    legend('Balance', 'PBAL', 'location', 'northeast') 
    legend('boxoff')
    title(filename)
   
  figure(3)
    plot(TIME, BPCX./BPIT, 'k', 'linewidth', 2, TIME, (BPLIM + BPSHI)./BPIT, 'r', 'linewidth', 2, TIME, BPTP./BPIT, 'b', 'linewidth', 2) 
    axis([TIME(1) TIME(end)])
    xlabel('time (s)')
    ylabel('Fraction of NBI power to')
    legend('CX', 'Orbot Loss + ShineThrou', 'To el/ion', 'location', 'northeast') 
    legend('boxoff')
    title(filename)
    %} 
endif




















	
