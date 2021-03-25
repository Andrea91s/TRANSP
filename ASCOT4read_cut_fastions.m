%%%%%%%%%%%%%%%%%5 HERE I READ SOME VARIABLE FROM TRANSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fast = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P17/29909P17_fi_1.cdf';
fast2 = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_fi_1.cdf';
state = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P17/29909P17_ps_ts1_state.cdf';
big = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P17/29909P17.CDF';
load('RBW.mat')
%CDF
TIME = ncread(big, 'TIME');
TIME3 = ncread(big, 'TIME3');
Xrho = ncread(big, 'X');
#OMEGA = ncread(big, 'OMEGA');
OMEGA = zeros(60,58); 
ND = ncread(big, 'ND');
NE = ncread(big, 'NE');
TI = ncread(big, 'TI');
TE = ncread(big, 'TE');
BDENS = ncread(big, 'BDENS');
DVOL = nc_read(big, 'DVOL');   % ZONE VOLUME
NTOT = ncread(fast2, 'NTOT_D_NBI')
fast = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P05/29909P05_fi_1.cdf';
RSURF = ncread(fast, 'RSURF');
ZSURF = ncread(fast, 'ZSURF');
%FAST IONS
%FAST IONS
NE_D_NBI = ncread(fast, 'NE_D_NBI');
NA_D_NBI = ncread(fast, 'NA_D_NBI');
E_D_NBI = ncread(fast, 'E_D_NBI');
A_D_NBI = ncread(fast, 'A_D_NBI');
NTHSURF = ncread(fast, 'NTHSURF');
R2D = ncread(fast, 'R2D');
Z2D = ncread(fast, 'Z2D');
X2D = ncread(fast, 'X2D');
BMVOL = ncread(fast, 'BMVOL');
BMVOLFLR = ncread(fast2, 'BMVOL');


F_D_NBI = ncread(fast, 'F_D_NBI');
F_D_NBIFLR = ncread(fast2, 'F_D_NBI');
NTOT_D_NBI = ncread(fast, 'NTOT_D_NBI');
NTOT_D_NBIFLR = ncread(fast2, 'NTOT_D_NBI');
TIMEfi = ncread(fast, 'TIME');
DT_AVG = ncread(fast, 'DT_AVG');
nzones = ncread(fast, 'nzones');
N_2DZONES = ncread(fast,'N_2DZONES');
NXSURF = ncread(fast,'NXSURF');
nsjdotb = ncread(fast,'nsjdotb');
nsnccw = ncread(fast,'nsnccw');
addpath('/home/andrea/Documents/MASTOrbit');
[EQD] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/g029909.00216',0);
%[EQD] = read_eqdsk('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_matlab/g29210.00255',0);
load('RBW.mat')

R = linspace(0.21, 1.43, 62);
Z = linspace(-1.01, 1.21, 112);

for i=1: length(R)-1
  Rhisto(i) = (R(i) + R(i+1))/2;
end
Rhisto = Rhisto.*100;

for i=1: length(Z)-1
  Zhisto(i) = (Z(i) + Z(i+1))/2;
end
Zhisto = Zhisto.*100;

PPP = linspace(0,1,121);
[Rascot, Zascot] = meshgrid(Rhisto, Zhisto);
R2D_Ao = reshape(Rascot,[],1);
Z2D_Ao = reshape(Zascot,[],1);
[inside,outside] = inpolygon(R2D_Ao,Z2D_Ao,RSURF(:,109), ZSURF(:,109));
[inside1,outside1] = inpolygon(R2D_Ao,Z2D_Ao,RSURF(:,97), ZSURF(:,97));
in = inside & ! outside;
out = inside1 & ! outside1;
out = ~out;
R2D_A = R2D_Ao(in & out);
Z2D_A = Z2D_Ao(in & out);

return

figure(1000)
scatter(R2D_Ao, Z2D_Ao, 'r')
hold all;
scatter(R2D_A,Z2D_A,'k')

close all;

BMVOL_Ao = (2*pi*( max(diff(Z2D_Ao))* max(diff(R2D_Ao)))*R2D_Ao);
BMVOL_A = BMVOL_Ao(in & out);

X2D_Ao = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2D_Ao, Z2D_Ao);
X2D_A = X2D_Ao(in & out);

ascot4GC = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A413.h5');
ascot4FO = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A404.h5');

dist4AGO = squeeze(ascot4GC.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adressogc = reshape(permute(dist4AGO, [4,3,2,1]), [75, 50, 6771]);
dist4adressgc = dist4adressogc(:,:,in & out)/4744;


dist4afo = squeeze(ascot4FO.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adressofo = reshape(permute(dist4afo, [4,3,2,1]), [75, 50, 6771]);
dist4adressfo = dist4adressofo(:,:,in & out)/5532;




pitch4 = ascot4FO.distributions.rzPitchEdist.abscissae.dim3;
energy4 = ascot4FO.distributions.rzPitchEdist.abscissae.dim4.*6.242e18;
for i=1: length(pitch4)-1
  p4(i) = (pitch4(i) + pitch4(i+1))/2;
endfor

for i=1: length(energy4)-1
  e4(i) = (energy4(i) + energy4(i+1))/2;
endfor





FIDA4GLOBALGC = zeros(75,50);
FIDA4GLOBALFO = zeros(75,50);



for i=1:length(BMVOL_A)
  b = dist4adressgc(:,:,i).*BMVOL_A(i);
  t = dist4adressfo(:,:,i).*BMVOL_A(i);
  FIDA4GLOBALGC = FIDA4GLOBALGC+b;
  FIDA4GLOBALFO = FIDA4GLOBALFO+t;
end 
FIDA4PGC = sum(FIDA4GLOBALGC,1).*(e4(2)-e4(1));
FIDA4EGC = sum(FIDA4GLOBALGC,2).*(p4(2)-p4(1));
FIDA4PFO = sum(FIDA4GLOBALFO,1).*(e4(2)-e4(1));
FIDA4EFO = sum(FIDA4GLOBALFO,2).*(p4(2)-p4(1));


figure(200)
subplot(1,3,1)
pcolor(p4,e4, FIDA4GLOBALGC)
shading interp
title('ASCOT4 GC, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading flat
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
%caxis([0 2e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)

subplot(1,3,2)
pcolor(p4,e4, FIDA4GLOBALFO)
shading interp
title('ASCOT4 FO, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading flat
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
%caxis([0 2e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)

subplot(1,3,3)
pcolor(p4,e4, 100.*(1-(FIDA4GLOBALFO./sum(FIDA4GLOBALFO(:)))./(FIDA4GLOBALGC./(sum(FIDA4GLOBALGC(:))))))
shading interp
title('Relative difference', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading flat
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
caxis([-100 100]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)



ASCOT4GC = sum((BMVOL_A.*squeeze(sum(squeeze(sum(dist4adressgc,1)),1))')).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2
ASCOT4FO = sum((BMVOL_A.*squeeze(sum(squeeze(sum(dist4adressfo,1)),1))')).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2

FIDA4PFO = sum(FIDA4GLOBALFO,1).*(e4(2)-e4(1));
FIDA4EFO = sum(FIDA4GLOBALFO,2).*(p4(2)-p4(1));
FIDA4PGC = sum(FIDA4GLOBALGC,1).*(e4(2)-e4(1));
FIDA4EGC = sum(FIDA4GLOBALGC,2).*(p4(2)-p4(1));

figure(300)
subplot(2,2,1)
plot(p4, FIDA4PFO,'r','linewidth',2,...
p4, FIDA4PGC,'b','linewidth',2)
legend('FO','GC')
xlabel('pitch')
ylabel('dFID(P)/dP density')


subplot(2,2,2)
plot(e4, FIDA4EFO,'r','linewidth',2,...
e4, FIDA4EGC,'b','linewidth',2)
legend('FO','GC')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')

subplot(2,2,3)
plot(p4, FIDA4PFO./sum(FIDA4PFO),'r','linewidth',2,...
p4, FIDA4PGC./sum(FIDA4PGC),'b','linewidth',2)
legend('FO','GC')
xlabel('pitch')
ylabel('dFID(P)/dP density')

subplot(2,2,4)
plot(e4, FIDA4EFO./sum(FIDA4EFO),'r','linewidth',2,...
e4,FIDA4EGC./sum(FIDA4EGC),'b','linewidth',2)
legend('FO','GC')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')


return
nccreate('AGO_fi_1.cdf','NE_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','NA_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','E_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4)});
nccreate('AGO_fi_1.cdf','A_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(p4)});
nccreate('AGO_fi_1.cdf','NTHSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','R2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('AGO_fi_1.cdf','Z2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('AGO_fi_1.cdf','X2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('AGO_fi_1.cdf','BMVOL', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('AGO_fi_1.cdf','RSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00121',121});
nccreate('AGO_fi_1.cdf','ZSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00121',121});
nccreate('AGO_fi_1.cdf','F_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4),'dim_00050',length(p4),'dim_04481',length(BMVOL_A)});
nccreate('AGO_fi_1.cdf','NTOT_D_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('AGO_fi_1.cdf','TIME', 'Format', '64bit', 'Datatype', 'double');
nccreate('AGO_fi_1.cdf','DT_AVG', 'Format', '64bit', 'Datatype', 'double');
nccreate('AGO_fi_1.cdf','nzones', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','N_2DZONES', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','NXSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','nsjdotb', 'Format', '64bit', 'Datatype', 'int32');
nccreate('AGO_fi_1.cdf','nsnccw', 'Format', '64bit', 'Datatype', 'int32');

ncwrite('AGO_fi_1.cdf','NE_D_NBI', length(e4));
ncwrite('AGO_fi_1.cdf','NA_D_NBI', length(p4));
ncwrite('AGO_fi_1.cdf','E_D_NBI', e4');
ncwrite('AGO_fi_1.cdf','A_D_NBI', p4');
ncwrite('AGO_fi_1.cdf','NTHSURF', NTHSURF);
ncwrite('AGO_fi_1.cdf','R2D', R2D_A);
ncwrite('AGO_fi_1.cdf','Z2D', Z2D_A);
ncwrite('AGO_fi_1.cdf','X2D', X2D_A);
ncwrite('AGO_fi_1.cdf','BMVOL', BMVOL_A);
ncwrite('AGO_fi_1.cdf','RSURF', RSURF);
ncwrite('AGO_fi_1.cdf','ZSURF', ZSURF);
ncwrite('AGO_fi_1.cdf','F_D_NBI', dist4adressfo);
ncwrite('AGO_fi_1.cdf','NTOT_D_NBI', ASCOT4FO);
ncwrite('AGO_fi_1.cdf','TIME', TIMEfi);
ncwrite('AGO_fi_1.cdf','DT_AVG',  DT_AVG);
ncwrite('AGO_fi_1.cdf','nzones', 60);
ncwrite('AGO_fi_1.cdf','N_2DZONES', length(BMVOL_A));
ncwrite('AGO_fi_1.cdf','NXSURF', NXSURF);
ncwrite('AGO_fi_1.cdf','nsjdotb', -1);
ncwrite('AGO_fi_1.cdf','nsnccw', -1);

