function [TIME, X, ZEFFP, ZEFFI] = nc_zeff_profile(filename)
% Script to access neutrons data from the special TRANSP output file
%
% INPUT
% filename	string		name of the special TRANSP output
%
% OUTPUT
% Zeff		array		radial coordinates of the non averaged flux surfaces in a given time in cm

% Read the normalized flux surface coordiinate
X = nc_read(filename, 'X');

% Read the time
TIME = nc_read(filename, 'TIME');   % time in sec

% Read the Zeff profile
ZEFFP = nc_read(filename, 'ZEFFP');

% Read the Zeff profile
ZEFFI = nc_read(filename, 'ZEFFI');