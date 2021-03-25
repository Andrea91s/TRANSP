function [FBEMIN FBEMAX FVPVMN FVPVMX FSHPER FSHWID TFSHON TFSHOF] = nc_fishbones_parameters(filename)
% Script to access neutrons data from TRANSP output file
%
% INPUT
% filename		name of the TRANSP file containing the results
%
% OUTPUT

% Read the total neutron emissivities
FBEMIN = nc_read(filename, 'FBEMIN');
FBEMAX = nc_read(filename, 'FBEMAX');
FVPVMN = nc_read(filename, 'FVPVMN');
FVPVMX = nc_read(filename, 'FVPVMX');
FBLTIM = nc_read(filename, 'FBLTIM');
FSHPER = nc_read(filename, 'FSHPER');
FSHWID = nc_read(filename, 'FSHWID');
TFSHON = nc_read(filename, 'TFSHON');
TFSHOF = nc_read(filename, 'TFSHOF');