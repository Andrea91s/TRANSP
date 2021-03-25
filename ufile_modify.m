clear
close all
UFILE = ufile_read(filename = 'RUNs/29880/U11/F29880.NERA', plotyn = 1, tslice = 0.22);
UFILEmod = UFILE;

% Set the density below to a certain value to that value
density_threshold = 1E12;
idx = find(UFILE.data <= density_threshold);
UFILEmod.data(idx) = density_threshold;


ufile_write(UFILEmod, 'RUNs/29880/U12/F29880.NERA');
UFILEmod = ufile_read(filename = 'RUNs/29880/U12/F29880.NERA', plotyn = 1, tslice = 0.22);





UFILE = ufile_read(filename = 'RUNs/29880/U11/F29880.ZFRA', plotyn = 1, tslice = 0.22);
UFILEmod = UFILE;

% Set Zeff flat
UFILEmod.data = 1.5*ones(size(UFILEmod.data));


ufile_write(UFILEmod, 'RUNs/29880/U12/F29880.ZFRA');
UFILEmod = ufile_read(filename = 'RUNs/29880/U12/F29880.ZFRA', plotyn = 1, tslice = 0.22);