% Modify Te to increase Yn to match MNEUT


filename = 'RUNs/26887/U01/F26887.TERA';

[UFILE] = ufile_read(filename, plotyn = 1, tslice = 0.25);
UFILE_MOD = UFILE;

% Set Zeff = min(Zeff) for all r/a < r/a for Zeff min 
UFILE_MOD.data = 1.05*UFILE_MOD.data;

% Write modified UFILE
filename_mod = 'RUNs/26887/U01/F26887.TERA_MOD';
ufile_write(UFILE_MOD, filename_mod);


% Read back and plot
ufile_read(filename_mod, plotyn = 1, tslice = 0.25);
