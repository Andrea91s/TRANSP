% Modify Zeff to increase Yn to match MNEUT


filename = 'RUNs/26887/U01/F26887.ZFRA';

[UFILE] = ufile_read(filename, plotyn = 1, tslice = 0.25);
UFILE_MOD = UFILE;

% Find the minimum in Zeff
for m=1:UFILE.M; 
  [Zeffmin(m), min_index(m)] = min(UFILE.data(m,:)); 
end




% Set Zeff = min(Zeff) for all r/a < r/a for Zeff min 
for m = 1:UFILE.M
  UFILE_MOD.data(m,1:min_index(m)) = Zeffmin(m);
end


% Write modified UFILE
filename_mod = 'RUNs/26887/U01/F26887.ZFRA_MOD';
ufile_write(UFILE_MOD, filename_mod);


% Read back and plot
ufile_read(filename_mod, plotyn = 1, tslice = 0.25);
