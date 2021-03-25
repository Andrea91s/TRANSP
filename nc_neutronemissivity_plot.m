function [rr, zz, ems, ynt] = nc_neutronemissivity_plot(TIME, R, Z, EMS, DV, time, fs)
% Function that return the neutron emissivity for a given time and plots it
%
% INPUT
% TIME          time array [s]
% R             time x zones R x poloidal angle [cm]
% Z             time x zones Z x poloidal angle [cm]
% EMS			neutron emissivity [cm-3 s-1]
% DV            time x zones volumes [cm^3]
% time          selected time for the plot [s]


% OUTPUT
% rr			R at selected time [cm]
% zz			Z at selected time [cm]
% ems			neutron emissivity [cm-3 s-1] at selected time
% ynt           neutron yield in [s-1]


% Find the index corresponding to the selected time
idx = max(find(TIME <= time));

% Select the time
rr = squeeze(R(idx,:,:));
zz = squeeze(Z(idx,:,:));
ems = squeeze(EMS(idx,:,:));

% volume
dv = DV(idx,:);

% Evaluate the total neutron yield
w = ems(:,1).*dv';
ynt = sum(w);

% Plot the results
%load MASTcmap2.dat
figure
pcolor(rr, zz, log10(ems'))
shading flat
if (fs == 1)
    hold on
	for k = 1:5:60
        plot(rr(k,:), zz(k,:), 'k')
    end
	hold off
endif

%set(gcf,'Colormap',MASTcmap2) 
xlabel('R (cm)', 'fontsize', 12)
ylabel('Z (cm)', 'fontsize', 12)
axis equal
axis([0 200 -150 150])
caxis([0 8])
h = colorbar;
set(h, 'title', 'log_{10} [cm^{-3} s^{-1}]')
title (['Neutron emissivity at t = ' num2str(TIME(idx)) ' s.'], 'fontsize', 14)
text(50, 130, sprintf('Neutron Yield = %1.4g s^{-1}', ynt))

