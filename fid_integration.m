function [fid, fin] = fid_integration(PA, EN, FID, BMVOL, id, PI1, PI2, EI1, EI2)
% Script that calculates the fast ions density in a point on the poloidal cross-section
% and in a user specifed energy region and range
%
% INPUT
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% BMVOL		array		zone volumes
% EI1		float		lower integral limit in energy
% EI2		float		upper integral limit in energy
% PI1		float		lower integral limit in pitch angle
% PI2		float		upper integral limit in pith angle
% id		integer		array index of FID which is closest to the requested point in the poloidal plane
%
% OUTPUT
% fid       float       fast ion density in units of (1/cm3) (1/4pi)
% fin       float       number of fast ion in units of (1/4pi)

% Example
% 


% Add the require path containing the routines
addpath('~/Documents/Octave');

% Fast ion distribution for the selected point in space
F = squeeze(FID(id,:,:));   % in units of (1/cm3) (1/eV) (1/4pi)

% Volume of the cell for the selected point in space
V = BMVOL(id);  % in units of cm3

% Note that due to the way the hintegral.m script is written
% the transpose of PA, EN and F should be passed to it!!
hi = hintegral2d(PA', EN', F', PI1, PI2, EI1, EI2);

% So the ouputs are
fid = hi;       % fast ion density in units of (1/cm3) (1/4pi)
fin = fid*V;     % number of fast ion in units of (1/4pi)

% Output the results
printf('For PA in [%f, %f] and EN in [%f, %f] keV:\n', PI1, PI2, EI1/1000, EI2/1000);
printf('Fast ion density: %1.4g [1/cm3 1/4pi]\n', fid);
printf('Fast ions: %1.4g [1/4pi]\n', fin);


rmpath('~/Documents/Octave');
