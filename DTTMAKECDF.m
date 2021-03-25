%%%%%%%%%%%%%%%%%5 HERE I READ SOME VARIABLE FROM TRANSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fasttbeam = '/media/andrea/My Passport Ultra/TRANSP/TRANSP_analysis_script/DTT/29270U24_fi_1.cdf';
fastfusn = '/media/andrea/My Passport Ultra/TRANSP/TRANSP_analysis_script/DTT/29270U22_fi_1.cdf';
state = '/media/andrea/My Passport Ultra/TRANSP/TRANSP_analysis_script/DTT/29270U24_ps_ts1_state.cdf';
big = '/media/andrea/My Passport Ultra/TRANSP/TRANSP_analysis_script/DTT/29270U24.CDF';

%CDF
TIME = ncread(big, 'TIME');
TIME3 = ncread(big, 'TIME3');
X = ncread(big, 'X');
OMEGA = ncread(big, 'OMEGA');
#ND = ncread(big, 'ND');
#NT = ncread(big, 'NT');
#TI = ncread(big, 'TI');
ND = 1e14*ones(60,376);
NT = 1e14*ones(60,376);
k=50;
TI = k*1000*ones(60,376);
return

%PLASMA STATE
R_grid = ncread(state, 'R_grid');
Z_grid = ncread(state, 'Z_grid');
BRRZ = ncread(state, 'BRRZ');
BphiRZ = ncread(state, 'BphiRZ');
BZRZ = ncread(state, 'BZRZ');


nccreate('TEST50keV.CDF','TIME', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME', length(TIME)});
nccreate('TEST50keV.CDF','TIME3', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME3', length(TIME3)});
nccreate('TEST50keV.CDF','X', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('TEST50keV.CDF','OMEGA', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('TEST50keV.CDF','ND', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('TEST50keV.CDF','TI', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('TEST50keV.CDF','NT', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});

nccreate('TEST50keV.CDF','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', 105});
nccreate('TEST50keV.CDF','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', 165});
nccreate('TEST50keV.CDF','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('TEST50keV.CDF','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('TEST50keV.CDF','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});

ncwrite('TEST50keV.CDF','TIME', TIME);
ncwrite('TEST50keV.CDF','TIME3', TIME3);
ncwrite('TEST50keV.CDF','X', X);
ncwrite('TEST50keV.CDF','OMEGA', OMEGA);
ncwrite('TEST50keV.CDF','ND', ND);
ncwrite('TEST50keV.CDF','NT', NT);
ncwrite('TEST50keV.CDF','TI', TI);


ncwrite('TEST50keV.CDF','R_grid', R_grid);
ncwrite('TEST50keV.CDF','Z_grid', Z_grid);
ncwrite('TEST50keV.CDF','BRRZ', BRRZ);
ncwrite('TEST50keV.CDF','BphiRZ', BphiRZ);
ncwrite('TEST50keV.CDF','BZRZ', BZRZ);
return

%FAST IONS
NE_D_NBI = ncread(fastfusn, 'NE_D_NBI');
NA_D_NBI = ncread(fastfusn, 'NA_D_NBI');
E_D_NBI = ncread(fastfusn, 'E_D_NBI');
A_D_NBI = ncread(fastfusn, 'A_D_NBI');
F_D_NBI = ncread(fastfusn, 'F_D_NBI');
NTOT_D_NBI = ncread(fastfusn, 'NTOT_D_NBI');
NE_D_NBI = ncread(fastfusn, 'NE_D_NBI');

NE_T_FUSN = ncread(fastfusn, 'NE_T_FUSN');
NA_T_FUSN = ncread(fastfusn, 'NA_T_FUSN');
E_T_FUSN = ncread(fastfusn, 'E_T_FUSN');
A_T_FUSN = ncread(fastfusn, 'A_T_FUSN');
F_T_FUSN = ncread(fastfusn, 'F_T_FUSN');
NTOT_T_FUSN = ncread(fastfusn, 'NTOT_T_FUSN');
NE_T_FUSN = ncread(fastfusn, 'NE_T_FUSN');

NE_T_NBI = ncread(fasttbeam, 'NE_T_NBI');
NA_T_NBI = ncread(fasttbeam, 'NA_T_NBI');
E_T_NBI = ncread(fasttbeam, 'E_T_NBI');
A_T_NBI = ncread(fasttbeam, 'A_T_NBI');
F_T_NBI = ncread(fasttbeam, 'F_T_NBI');
NTOT_T_NBI = ncread(fasttbeam, 'NTOT_T_NBI');
NE_T_NBI = ncread(fasttbeam, 'NE_T_NBI');

NTHSURF = ncread(fastfusn, 'NTHSURF');
R2D = ncread(fastfusn, 'R2D');
Z2D = ncread(fastfusn, 'Z2D');
X2D = ncread(fastfusn, 'X2D');
BMVOL = ncread(fastfusn, 'BMVOL');
RSURF = ncread(fastfusn, 'RSURF');
ZSURF = ncread(fastfusn, 'ZSURF');
F_D_NBI = ncread(fastfusn, 'F_D_NBI');
NTOT_D_NBI = ncread(fastfusn, 'NTOT_D_NBI');
TIMEfi = ncread(fastfusn, 'TIME');
DT_AVG = ncread(fastfusn, 'DT_AVG');
nzones = ncread(fastfusn, 'nzones');
N_2DZONES = ncread(fastfusn,'N_2DZONES');
NXSURF = ncread(fastfusn,'NXSURF');
nsjdotb = ncread(fastfusn,'nsjdotb');
nsnccw = ncread(fastfusn,'nsnccw');


return




nccreate('TEST_fi.cdf','NE_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NA_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NE_T_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NA_T_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NE_T_FUSN', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NA_T_FUSN', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','E_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_D_NBI)});
nccreate('TEST_fi.cdf','A_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(A_D_NBI)});
nccreate('TEST_fi.cdf','E_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_T_NBI)});
nccreate('TEST_fi.cdf','A_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(A_T_NBI)});
nccreate('TEST_fi.cdf','E_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_T_FUSN)});
nccreate('TEST_fi.cdf','A_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(A_T_FUSN)});

nccreate('TEST_fi.cdf','F_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_D_NBI),'dim_00050',length(A_D_NBI),'dim_00840',length(BMVOL)});
nccreate('TEST_fi.cdf','F_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_T_NBI),'dim_00050',length(A_T_NBI),'dim_00840',length(BMVOL)});
nccreate('TEST_fi.cdf','F_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00100',length(E_T_FUSN),'dim_00050',length(A_T_FUSN),'dim_00840',length(BMVOL)});

nccreate('TEST_fi.cdf','NTOT_D_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('TEST_fi.cdf','NTOT_T_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('TEST_fi.cdf','NTOT_T_FUSN', 'Format', '64bit', 'Datatype', 'double');

nccreate('TEST_fi.cdf','NTHSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','R2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00840', length(BMVOL)});
nccreate('TEST_fi.cdf','Z2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00840', length(BMVOL)});
nccreate('TEST_fi.cdf','X2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_060840', length(BMVOL)});
nccreate('TEST_fi.cdf','BMVOL', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_06771', length(BMVOL)});
nccreate('TEST_fi.cdf','RSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00041',41});
nccreate('TEST_fi.cdf','ZSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00041',41});
nccreate('TEST_fi.cdf','TIME', 'Format', '64bit', 'Datatype', 'double');
nccreate('TEST_fi.cdf','DT_AVG', 'Format', '64bit', 'Datatype', 'double');
nccreate('TEST_fi.cdf','nzones', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','N_2DZONES', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','NXSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','nsjdotb', 'Format', '64bit', 'Datatype', 'int32');
nccreate('TEST_fi.cdf','nsnccw', 'Format', '64bit', 'Datatype', 'int32');



ncwrite('TEST_fi.cdf','NE_D_NBI', NE_D_NBI);
ncwrite('TEST_fi.cdf','NE_T_NBI', NE_T_NBI);
ncwrite('TEST_fi.cdf','NE_T_FUSN', NE_T_FUSN);

ncwrite('TEST_fi.cdf','NA_D_NBI', NA_D_NBI);
ncwrite('TEST_fi.cdf','NA_T_NBI', NA_T_NBI);
ncwrite('TEST_fi.cdf','NA_T_FUSN', NA_T_FUSN);

ncwrite('TEST_fi.cdf','E_D_NBI', E_D_NBI);
ncwrite('TEST_fi.cdf','E_T_NBI', E_T_NBI);
ncwrite('TEST_fi.cdf','E_T_FUSN', E_T_FUSN);

ncwrite('TEST_fi.cdf','A_D_NBI', A_D_NBI);
ncwrite('TEST_fi.cdf','A_T_NBI', A_T_NBI);
ncwrite('TEST_fi.cdf','A_T_FUSN', A_T_FUSN);

ncwrite('TEST_fi.cdf','F_D_NBI', F_D_NBI);
ncwrite('TEST_fi.cdf','F_T_NBI', F_T_NBI);
ncwrite('TEST_fi.cdf','F_T_FUSN', F_T_FUSN);

ncwrite('TEST_fi.cdf','NTOT_D_NBI', NTOT_D_NBI);
ncwrite('TEST_fi.cdf','NTOT_T_NBI', NTOT_T_NBI);
ncwrite('TEST_fi.cdf','NTOT_T_FUSN', NTOT_T_FUSN);

ncwrite('TEST_fi.cdf','NTHSURF', NTHSURF);
ncwrite('TEST_fi.cdf','R2D', R2D);
ncwrite('TEST_fi.cdf','Z2D', Z2D);
ncwrite('TEST_fi.cdf','X2D', X2D);
ncwrite('TEST_fi.cdf','BMVOL', BMVOL);
ncwrite('TEST_fi.cdf','RSURF', RSURF);
ncwrite('TEST_fi.cdf','ZSURF', ZSURF);


ncwrite('TEST_fi.cdf','TIME', TIMEfi);
ncwrite('TEST_fi.cdf','DT_AVG',  DT_AVG);
ncwrite('TEST_fi.cdf','nzones', 60);
ncwrite('TEST_fi.cdf','N_2DZONES', length(BMVOL));
ncwrite('TEST_fi.cdf','NXSURF', NXSURF);
ncwrite('TEST_fi.cdf','nsjdotb', 1);
ncwrite('TEST_fi.cdf','nsnccw', 1);



 