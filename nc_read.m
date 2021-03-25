function [data] = nc_read(filename, variable)
% Function that read a variable in a NetCDF file
%
% INPUT
%
% filename		NetCDF filename
% variable		variable to read

% Read the NETCDF file
nc = netcdf(filename,'r','nowrite');
data = nc{variable}(:);

