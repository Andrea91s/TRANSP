function [r, rems] = nc_neutronemissivity_onaxis(TIME, R, Z, EMS, RAXIS, ZAXIS, time, plotyn)
% Function that interpolates the neutron emissivity along the major radius
%
% INPUTS
%
% TIME          time array [s]
% R             radial coordinates [cm]
% Z             vertical coordinates [cm]
% EMS           neutron emissivity in [cm^-3 s^-1]
% RAXIS         magnetic axis R [cm]
% ZAXIS         magnetic axis Z [cm]
% time          time for which the 1D radial profile is required [s]
% plotyn        keyword
%
% OUTPUTS
% r             radial coordinate at t = time for Z = ZAXIS(time)
% rems          radial profile of the neutron emissivity [cm^-3 s^-1]

% Define the radial coordinate
r = linspace(0, 200, 1000);

% Find the index of the corresponding time
idx = max(find(TIME <= time));

% Extract the data at the selected time
rr = squeeze(R(idx,:,:)); 
zz = squeeze(Z(idx,:,:)); 
ems = squeeze(EMS(idx,:,:));

% Calculates the interpolation
rems = griddata(rr, zz, ems, r, ZAXIS(idx));
rems(find(isnan(rems))) = 0;

% Makes the plot
if (plotyn == 1)
    figure
    plot(r, rems, 'linewidth', 2)
    xlabel('major radius (cm)', 'fontsize', 14)
    ylabel('Neutron Yield (cm^{-3} s^{-1})', 'fontsize', 14)
    title(['Profile at Z = ' num2str(ZAXIS(idx)) ' cm, for t = ' num2str(TIME(idx)) ' s.'])
endif
