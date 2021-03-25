% Calculations of the expected count rates at the neutron camera
% for given MAST U Scenario and RNC design
close all;
clear all;

% Load the MAST Upgrade neutron emissivity profile from TRANSP
% predictive runs (the scenario is specified within the script
% called below). 
TRANSP_Analysis;
return
clear -x YN;


% Generate the RNC geometry
addpath('/home/andrea/OFF_AXIS_NCU/');
    %NCU_LOS_06CH;
    %NCU_LOS_06CH_3_3_deg;
    NCU_OFF;
    %JET_LOS_19CH_out;

    clear -x YN RNC TNC NC
rmpath('/home/andrea/OFF_AXIS_NCU/');


% Calculates the l
addpath('/home/andrea//NeutronCameraUpgrade/LIN/');
    [LINE] = line_integrated_emissivity(EQD = 1, YN, RNC, plotyn = 0);
    [TLINE] = line_integrated_emissivity(EQD = 1, YN, TNC, plotyn = 0);
rmpath('/home/andrea/NeutronCameraUpgrade/LIN/');

% Plot the counts for LINE and TLINE showing the expected profile
% for two repeated pulses in which the camera has been rotated by 
% 1 deg
figure(10);
plot(RNC.ImpactParameter.p, LINE.Counts, 'bo', 'markersize', 10, 'linewidth', 2, TNC.ImpactParameter.p, TLINE.Counts, 'ro', 'markersize', 10, 'linewidth', 2)
axis([0.4 1.2 0 1.2*max(LINE.Counts)])
legend('Ref, position', '1 DEG rotated')
h = xlabel('Impact parameter (m)'); set(h , 'fontsize', 12)
h = ylabel('Counts'); set(h , 'fontsize', 12)
h = title(strrep(YN.Scenario, '_', ' ')); set(h , 'fontsize', 12)

print -dpng '9999.png'