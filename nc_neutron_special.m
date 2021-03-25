function [RNA, ZNA, BBN4, BTN4] = nc_neutron_special(filename)
% Script to access neutrons data from the special TRANSP output file
%
% INPUT
% filename	string		name of the special TRANSP output
%
% OUTPUT
% RNA		array		radial coordinates of the non averaged flux surfaces in a given time in cm
% ZNA		array		vertical coordinates of the non averaged flux surfaces in a given time in cm
% BBN4		array		non-flux averaged beam-beam neutron emissivity in cm-3 s-1
% BTN4		array		non-flux averaged beam-thermal neutron emissivity in cm-3 s-1


% Read the non-flux surface averaged neutron emissivities for the BB and BT components
RNA = nc_read(filename, 'R');
ZNA = nc_read(filename, 'Z');
BBN4 = nc_read(filename, 'BBN4');
BTN4 = nc_read(filename, 'BTN4');

