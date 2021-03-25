
% Name of the file from which the data should be read
filename = '/home/andrea/Documents/TRANSP/29880/29880U16_fi_4.cdf';
%filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U16/29880U16_fi_2.cdf';
%filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U31/29880U31_fi_2.cdf';
%filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U32/29880U32_fi_2.cdf';
%filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U32/29880U32_fi_1.cdf';
%filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U31/29880U31_fi_1.cdf';

% Read the data for the selected point in space
%[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 90, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
%[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 100, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
%[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 110, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 120, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);

% Plot the boundaries
% fid_boundary(PA, EN, R2D, Z2D, FID, id, IE, CM);


% Evaluate the integrals of the FID in the different regions as
% indicated in the boundary plot above

% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s and 0.263 s, R = 0.9 m Z = 0.0 m
% ----------------------------------------------------------------------
%{
% Region 1
PI1 = -1;
PI2 = 0.05;
EI1 = 0;
EI2 = 6E4;
[fid1, fin1] = fid_integration(PA, EN, FID, BMVOL, id, PI1, PI2, EI1, EI2);
clear PI1 PI2 EI1 EI2
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2

% Region 2
PI1 = 0.05*ones(1, 5);
PI2 = [0.25 0.35 0.45 0.55 0.65];
EI1 = [0 7.5 12.5 32.5 57.5]*1E3;
EI2 = [7.5 12.5 32.5 57.5 60]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2


% Region 3
PI1 = [0.25 0.35 0.45 0.55 0.65];
PI2 = ones(1, 5);
EI1 = [0 7.5 12.5 32.5 57.5]*1E3;
EI2 = [7.5 12.5 32.5 57.5 60]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2


% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t

% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt
%}

%{
% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s and 0.263 s, R = 1.0 m Z = 0.0 m
% ----------------------------------------------------------------------

% The parameter DL below is used to change the boundary between trapped and
% passing fast ions by an arbitrary amount simulating the uncertainty in the
% boundaries. If  not used in this sense it should be set to zero! A typical
% value used for the error analysis is 0.1
DL = 0.0;
% Region 1 - Counter Passing
PI1 = [-1.00 -1.00 -0.85 -0.75];
PI2 = [-0.35 -0.45 -0.55 -0.65];
EI1 = [ 0.00  7.50 12.50 17.50]*1E3;
EI2 = [ 7.5  12.50 17.50 22.50]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region 2 - Stagnation
PI1 = [-1.00 -1.00 -1.00 -0.65 -0.55 -0.45 -0.35];
PI2 = [-0.85 -0.75  0.05  0.05  0.05  0.05  0.05];
EI1 = [12.50  17.5 22.50 17.50 12.50  7.50  0.00]*1E3;
EI2 = [17.50  22.5 60.00 22.50 17.50 12.50  7.50]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2


% Region 3 - Trapped
PI1 = [ 0.05  0.05  0.05  0.05  0.05];
PI2 = [ 0.55  0.65  0.55  0.65  0.75]+DL;
EI1 = [ 0.00 22.50 27.50 32.50 57.50]*1E3;
EI2 = [22.50 27.50 32.50 57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2


% Region 4 - Passing
PI1 = [ 0.55  0.65  0.55  0.65  0.75]+DL;
PI2 = [ 1.00  1.00  1.00  1.00  1.00];
EI1 = [ 0.00 22.50 27.50 32.50 57.50]*1E3;
EI2 = [22.50 27.50 32.50 57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2

% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t
fid4t

% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt
(fid3t+fid4t)/fidt
fid4t + fid3t
%}


%{
% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s, R = 1.0 m Z = 0.0 m
% with TF ripples on
% ----------------------------------------------------------------------

% The parameter DL below is used to change the boundary between trapped and
% passing fast ions by an arbitrary amount simulating the uncertainty in the
% boundaries. If  not used in this sense it should be set to zero! A typical
% value used for the error analysis is 0.1
DL = 0.0;
% Region 1 - Counter Passing
PI1 = [-1.00 -1.00 -1.00 -0.85];
PI2 = [-0.35 -0.45 -0.55 -0.65];
EI1 = [ 0.00  7.50 12.50 17.50]*1E3;
EI2 = [ 7.5  12.50 17.50 22.50]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region 2 - Stagnation
PI1 = [-1.00 -0.65 -1.00 -0.55 -0.45 -0.35];
PI2 = [-0.85  0.05  0.05  0.05  0.05  0.05];
EI1 = [17.50  17.5 22.50 12.50  7.50  0.00]*1E3;
EI2 = [22.50  22.5 60.00 17.50 12.50  7.50]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2


% Region 3 - Trapped
PI1 = [ 0.05  0.05  0.05  0.05];
PI2 = [ 0.45  0.55  0.65  0.75]+DL;
EI1 = [ 0.00  7.50 22.50 57.50]*1E3;
EI2 = [ 7.50 22.50 57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2


% Region 4 - Passing
PI1 = [ 0.45  0.55  0.65  0.75]+DL;
PI2 = [ 1.00  1.00  1.00  1.00];
EI1 = [ 0.00  7.50 22.50 57.50]*1E3;
EI2 = [ 7.50 22.50 57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2

% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t
fid4t

% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt
fid1t + fid2t + fid3t + fid4t
(fid3t+fid4t)/fidt
fid4t + fid3t

% and corresponding boundaries     
B1.x = [-1.05 -0.85 -0.85 -0.65 -0.65 -0.55 -0.55 -0.45 -0.45 -0.35 -0.35];
B1.y = [17.50 17.50 22.50 22.50 17.50 17.50 12.50 12.50  7.50  7.50  0.00]*1E3;
B2.x = [0.05 0.05];
B2.y = [0.00 62.5]*1E3;
B3.x = [0.45 0.45 0.55 0.55 0.65 0.65 0.75 0.75];
B3.y = [0.00 7.50 7.50 22.5 22.5 57.5 57.5 62.5]*1E3;
%}


%{
% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s, R = 1.1 m Z = 0.0 m
% ----------------------------------------------------------------------
% Region counter passin
PI1 = [-1 -1];
PI2 = [-0.45 -0.35];
EI1 = [0 12.5]*1E3;
EI2 = [12.5 60]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region stagnation
PI1 = [-0.35 -0.35 -0.15 -0.05];
PI2 = [-0.25 0.05 0.05 0.05];
EI1 = [17.5 22.5 7.5 0]*1E3;
EI2 = [22.5 60 22.5 7.5]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2

% Region trapped
PI1 = [-0.45 -0.35 -0.25 -0.15 0.05 0.65];
PI2 = [-0.15 -0.15 -0.15 -0.05 0.65 0.75];
EI1 = [0 12.5 17.5 0 0 22.5]*1E3;
EI2 = [12.5 17.5 22.5 7.5 60 60]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2

% Region co-passing
PI1 = [0.65 0.75];
PI2 = [1 1];
EI1 = [0 22.5]*1E3;
EI2 = [22.5 60]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2


% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);

% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt
%}

%{
% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.263 s, R = 1.1 m Z = 0.0 m
% ----------------------------------------------------------------------
% Region counter passin
PI1 = [-1 -1];
PI2 = [-0.45 -0.35];
EI1 = [0 7.5]*1E3;
EI2 = [7.5 60]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region stagnation
PI1 = [-0.35 -0.15 -0.05];
PI2 = [ 0.05  0.05  0.05];
EI1 = [17.50  7.50  0.00]*1E3;
EI2 = [60.00 17.50  7.50]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2

% Region trapped
PI1 = [-0.45 -0.35  0.05  0.65];
PI2 = [-0.05 -0.15  0.65  0.75];
EI1 = [ 0.00  7.50  0.00 22.50]*1E3;
EI2 = [ 7.50 17.50 60.00 60.00]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2

% Region co-passing
PI1 = [0.65 0.75];
PI2 = [1 1];
EI1 = [0 22.5]*1E3;
EI2 = [22.5 60]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2


% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t
fid4t


% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt
%}



%{
% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.255 and 0.263 s, R = 1.2 m Z = 0.0 m
% ----------------------------------------------------------------------
% Region counter passin
PI1 = [-1.00 -1.00 -1.00 -1.00];
PI2 = [-0.65 -0.55 -0.45 -0.35];
EI1 = [ 0.00  7.5  22.50 57.50]*1E3;
EI2 = [ 7.50 22.50 57.50 60.00]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region potato
PI1 = [-0.55 -0.45];
PI2 = [-0.45 -0.35];
EI1 = [17.50 42.50]*1E3;
EI2 = [22.50 52.50]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2


% Region trapped
PI1 = [-0.65 -0.55 -0.45 -0.35 -0.35 -0.45 -0.35  0.05  0.05  0.05  0.05];
PI2 = [-0.05 -0.05 -0.05 -0.05 -0.15 -0.15 -0.15  0.75  0.85  0.75  0.65];
EI1 = [ 0.00  7.50 17.50 42.50 47.50 52.50 57.50  0.00 27.50 32.50 52.50 ]*1E3;
EI2 = [ 7.50 17.50 42.50 47.50 52.50 57.50 62.50 27.50 32.50 52.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2

% Region stagnation
PI1 = [-0.05 -0.15];
PI2 = [ 0.05  0.05];
EI1 = [ 0.00 47.50]*1E3;
EI2 = [47.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2

% Region co-passing
PI1 = [0.75 0.85 0.95];
PI2 = [1 1 1];
EI1 = [0 27.50 57.5]*1E3;
EI2 = [27.5 57.5 62.5]*1E3;
for k = 1:length(PI1)
    [fid5(k), fin5(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid5t = sum(fid5);
fin5t = sum(fin5);
clear PI1 PI2 EI1 EI2

% Region lost
PI1 = [0.65 0.65 0.75];
PI2 = [0.95 0.85 0.85];
EI1 = [57.5 52.50 32.5]*1E3;
EI2 = [62.5 57.50  52.50]*1E3;
for k = 1:length(PI1)
    [fid6(k), fin6(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid6t = sum(fid6);
fin6t = sum(fin6);
clear PI1 PI2 EI1 EI2


% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t
fid4t
fid5t
fid6t


% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt
fid5t/fidt
fid6t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt + fid5t/fidt + fid6t/fidt
%}




% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.255 and 0.258 s, R = 1.2 m Z = 0.0 m
% with TF ripples
% ----------------------------------------------------------------------
% Region counter passin
PI1 = [-1.00 -1.00 -1.00 -1.00];
PI2 = [-0.65 -0.55 -0.45 -0.35];
EI1 = [ 0.00  7.50 27.50 57.50]*1E3;
EI2 = [ 7.50 27.50 57.50 60.00]*1E3;
for k = 1:length(PI1)
    [fid1(k), fin1(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid1t = sum(fid1);
fin1t = sum(fin1);
clear PI1 PI2 EI1 EI2


% Region potato
PI1 = [-0.55 -0.45];
PI2 = [-0.45 -0.35];
EI1 = [17.50 42.50]*1E3;
EI2 = [27.50 57.50]*1E3;
for k = 1:length(PI1)
    [fid2(k), fin2(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid2t = sum(fid2);
fin2t = sum(fin2);
clear PI1 PI2 EI1 EI2


% Region trapped
PI1 = [-0.65 -0.55 -0.45 -0.35 -0.35  0.05  0.05  0.05  0.05  0.05];
PI2 = [-0.05 -0.05 -0.05 -0.05 -0.15  0.75  0.85  0.75  0.65  0.55];
EI1 = [ 0.00  7.50 17.50 42.50 57.50  0.00 27.50 32.50 47.50 52.50]*1E3;
EI2 = [ 7.50 17.50 42.50 57.50 62.50 27.50 32.50 47.50 52.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid3(k), fin3(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid3t = sum(fid3);
fin3t = sum(fin3);
clear PI1 PI2 EI1 EI2

% Region stagnation
PI1 = [-0.05 -0.15];
PI2 = [ 0.05  0.05];
EI1 = [ 0.00 57.50]*1E3;
EI2 = [57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid4(k), fin4(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid4t = sum(fid4);
fin4t = sum(fin4);
clear PI1 PI2 EI1 EI2

% Region co-passing
PI1 = [ 0.75 0.85   0.95];
PI2 = [ 1.00 1.00   1.00];
EI1 = [ 0.00 27.50 57.50]*1E3;
EI2 = [27.50 57.50 62.50]*1E3;
for k = 1:length(PI1)
    [fid5(k), fin5(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid5t = sum(fid5);
fin5t = sum(fin5);
clear PI1 PI2 EI1 EI2

% Region lost
PI1 = [ 0.55  0.55  0.65  0.75];
PI2 = [ 0.95  0.85  0.85  0.85];
EI1 = [57.50 52.50 47.50  32.5]*1E3;
EI2 = [62.50 57.50 52.50 47.50]*1E3;
for k = 1:length(PI1)
    [fid6(k), fin6(k)] = fid_integration(PA, EN, FID, BMVOL, id, PI1(k), PI2(k), EI1(k), EI2(k));
end
fid6t = sum(fid6);
fin6t = sum(fin6);
clear PI1 PI2 EI1 EI2


% Total integral
fidt = fid_integration(PA, EN, FID, BMVOL, id, -1, 1, 0, 6E4);
fid1t
fid2t
fid3t
fid4t
fid5t
fid6t


% Fractions
fid1t/fidt
fid2t/fidt
fid3t/fidt
fid4t/fidt
fid5t/fidt
fid6t/fidt

fid1t/fidt + fid2t/fidt + fid3t/fidt + fid4t/fidt + fid5t/fidt + fid6t/fidt





% ***************************************************************
% Figure
figure(1);
      pcolor(U.PA, U.EN, squeeze(U.FID(U.id,:,:))');
      shading flat
      xlabel('pitch angle')
      ylabel('energy eV')
      axis([-1 1 0 6E4])
      title(['FID at R = ' num2str(U.R2D(U.id)) ' cm, Z = ' num2str(U.Z2D(U.id)) ' cm.'])
      colorbar

      
return
      
B1.x = [-1.05 -0.85 -0.85 -0.75 -0.75 -0.65 -0.65 -0.55 -0.45 -0.45 -0.35 -0.35];
B1.y = [12.50 12.50 17.50 17.50 22.50 22.50 17.50 12.50 12.50  7.50  7.50  0.00]*1E3;
B2.x = [0.05 0.05];
B2.y = [0.00 62.5]*1E3;
B3.x = [0.55 0.55 0.65 0.65 0.55 0.55 0.65 0.65 0.75 0.75]+DL;
B3.y = [0.00 22.5 22.5 27.5 27.5 32.5 32.5 57.5 57.5 62.5]*1E3;     

% -------------------------------------------------------
% Orbit Topology for 29880 at 0.255 s, R = 1.0, Z = 0.0 m
% -------------------------------------------------------
B1.x = [-1.05 -0.85 -0.85 -0.75 -0.75 -0.65 -0.65 -0.55 -0.45 -0.45 -0.35 -0.35];
B1.y = [12.50 12.50 17.50 17.50 22.50 22.50 17.50 12.50 12.50  7.50  7.50  0.00]*1E3;
B2.x = [0.05 0.05];
B2.y = [0.00 62.5]*1E3;
B3.x = [0.55 0.55 0.65 0.65 0.55 0.55 0.65 0.65 0.75 0.75];
B3.y = [0.00 22.5 22.5 27.5 27.5 32.5 32.5 57.5 57.5 62.5]*1E3;


% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.255 and 0.263 s, R = 1.2 m Z = 0.0 m
% ----------------------------------------------------------------------
B1.x = [-0.65 -0.65 -0.55 -0.55 -0.45 -0.45 -0.35 -0.35];
B1.y = [ 0.00  7.50  7.50 22.50 22.50 57.50 57.50 62.50]*1E3;
B2.x = [-0.55 -0.55 -0.45 -0.45 -0.55];
B2.y = [17.50 22.50 22.50 17.50 17.50]*1E3;
B3.x = [-0.45 -0.45 -0.35 -0.35 -0.45];
B3.y = [42.50 52.50 52.50 42.50 42.50]*1E3;
B4.x = [-0.15 -0.15 -0.05 -0.05];
B4.y = [62.50 47.50 47.50  0.00]*1E3;
B5.x = [0.05  0.05];
B5.y = [0.00 62.50]*1E3;
B6.x = [0.65 0.65 0.75 0.75 0.85 0.85 0.95 0.95];
B6.y = [62.50 52.50 52.50 32.50 32.50 57.50 57.50 62.50]*1E3;
B7.x = [0.75 0.75 0.85 0.85 0.95 0.95];
B7.y = [0.00 27.50 27.50 57.50 57.50 62.50]*1E3;

% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.255 and 0.263 s, R = 1.2 m Z = 0.0 m
% with TF ripples
% ----------------------------------------------------------------------
B1.x = [-0.65 -0.65 -0.55 -0.55 -0.45 -0.45 -0.35 -0.35];
B1.y = [ 0.00  7.50  7.50 27.50 27.50 57.50 57.50 62.50]*1E3;
B2.x = [-0.55 -0.55 -0.45 -0.45 -0.55];
B2.y = [17.50 27.50 27.50 17.50 17.50]*1E3;
B3.x = [-0.45 -0.45 -0.35 -0.35 -0.45];
B3.y = [42.50 57.50 57.50 42.50 42.50]*1E3;
B4.x = [-0.15 -0.15 -0.05 -0.05];
B4.y = [62.50 57.50 57.50  0.00]*1E3;
B5.x = [0.05  0.05];
B5.y = [0.00 62.50]*1E3;
B6.x = [0.55   0.55  0.65  0.65  0.75  0.75  0.85  0.85  0.95  0.95];
B6.y = [62.50 52.50 52.50 47.50 47.50 32.50 32.50 57.50 57.50 62.50]*1E3;
B7.x = [0.75  0.75  0.85  0.85  0.95  0.95];
B7.y = [0.00 27.50 27.50 57.50 57.50 62.50]*1E3;

hold on
h = stairs(B1.x, B1.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B2.x, B2.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B3.x, B3.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B4.x, B4.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B5.x, B5.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B6.x, B6.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
h = stairs(B7.x, B7.y);
set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 2)
hold off
