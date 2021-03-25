function [EMS, UEMS] = nc_neutronemissivity(filename, variable, NTHETA, method);
% Function that returns the neutron emissivity from the NetCDF file
% specified in the
%
% INPUT
% filename		NetCDF file name
% variable		specifies the neutron emissivity component
%				BBNTX_DD:long_name = "DD BEAM-BEAM NEUTRONS 
%				BTNTX_DD:long_name = "DD BEAM-TARGET NEUTRONS
%				BTNTX:long_name = "BEAM-TARGET NEUTRONS
%				BBNTX:long_name = "BEAM-BEAM NEUTRONS
%				THNTX:long_name = "THERMONUCLEAR NEUTRONS
%				TTNTX:long_name = "TOTAL NEUTRONS
%				FTOTDT:long_name = "TOTAL D-T FUSION
%				FTOTDDN:long_name = "TOTAL D(D,N)HE3 FUSION 
%				FTOTDDP:long_name = "TOTAL D(D,P)T FUSION
%				THNTX_DD:long_name = "DD THERMONUCLEAR NEUTRONS 
% NTHETA        number of poloidal angles used to create the 2D map
%
% Method:       0   Add a point on the magnetic axis to avoid the hole
%                   in the 2D pcolor plot of the neutron emissivity
%
%               1   Do not add the point on the magnetic axis 
%
% OUTPUT
% EMS			neutron emissivity time x length(R) x length(theta) coordinates in cm-3 s-1
% UEMS			neutron emissivity time x length(R) + 1 x length(theta) coordinates in cm-3 s-1


% Read the total neutron emissivity
data = nc_read(filename, variable); % total neutrons n/s/cm3

% Create an empty array
[m, n] = size(data);

% ********************************************************
% Routine that add one fake flux surface at the magnetic
% axis to avoid the void point in the centre
% ********************************************************
if (method == 0)
    n = n + 1;	% add one flux sruface for the magnetic axis
    %k = 50;     % number of polodal angles
    EMS = zeros(m, n, NTHETA);

    % the first flux surface contains the magnetic axis
    for i = 1:m
    EMS(i,1,:) = data(i,1);
    %EMS(i,1,:) = data_interp(i,1);
    end

    % Averaged neutron emissivity
    for i = 1: m
        for j = 2:n
            EMS(i,j,:) = data(i,j-1);
        end
    end
end

% ********************************************************
% Routine that does NOT one fake flux surface at the magnetic
% axis to avoid the void point in the centre
% ********************************************************
if (method == 1)
    EMS = zeros(m, n, NTHETA);

    % Averaged neutron emissivity
    for i = 1: m
        for j = 1:n
            EMS(i,j,:) = data(i,j);
        end
    end
end

UEMS = zeros(m, n+1, NTHETA);
UEMS(:,1,:) = EMS(:,1,:);
UEMS(:,2:n+1,:) = EMS(:,1:n,:);



