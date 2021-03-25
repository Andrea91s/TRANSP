#DEMO = load('/home/andrea/DEMO/ascot.h5');
#addpath('/home/andrea/Documents/MASTOrbit');
#[EQD] = read_eqdsk('/home/andrea/DEMO/Input_ASCOT/Equil_AR3d1_2015_04_v2_SOF_CSred_fine_final.eqdsk',0);
load('RBW.mat')
%%%%%%%%%%CREATE THE BIG CDF.FILE

TIME = linspace(1,10,10)';
TIME3 = linspace(1,10,10)';
ND = (DEMO.plasma._1d.ni(1,:)'./1e6).*ones(1,length(TIME));
NT = DEMO.plasma._1d.ni(2,:)'./1e6.*ones(1,length(TIME));
NE = DEMO.plasma._1d.ne'./1e6.*ones(1,length(TIME));
TI = DEMO.plasma._1d.ti'.*ones(1,length(TIME));
TE = DEMO.plasma._1d.te'.*ones(1,length(TIME));
X  = DEMO.plasma._1d.rho'.*ones(1,length(TIME));
OMEGA = ones(length(X),1).*ones(1,length(TIME));



R_grid = DEMO.bfield.r';
Z_grid = DEMO.bfield.z';
BRRZ = load('/home/andrea/DEMO/Input_ASCOT/DEMO_BR.dat')';
BphiRZ = load('/home/andrea/DEMO/Input_ASCOT/DEMO_BPHI.dat')';
BZRZ = load('/home/andrea/DEMO/Input_ASCOT/DEMO_BZ.dat')';





  %THIS IS THE HARDEST PART, WILL DO IT FOR LAST

NXSURF = 50;
RSURF = ones(length(EQD.LCFS_R), NXSURF).*EQD.LCFS_R.*100;
ZSURF = ones(length(EQD.LCFS_Z), NXSURF).*EQD.LCFS_Z.*100;
NTHSURF = length(EQD.LCFS_R);

R = DEMO.distributions.rzPitchEdist.abscissae.dim1;
Z = DEMO.distributions.rzPitchEdist.abscissae.dim2;
P = DEMO.distributions.rzPitchEdist.abscissae.dim3;
E = DEMO.distributions.rzPitchEdist.abscissae.dim4.*6.242e18;

for i=1: length(R)-1
  Rhisto(i) = (R(i) + R(i+1))/2;
end
Rhisto = Rhisto.*100;

for i=1: length(Z)-1
  Zhisto(i) = (Z(i) + Z(i+1))/2;
end
Zhisto = Zhisto.*100;

[Rascot, Zascot] = meshgrid(Rhisto, Zhisto);


R2Do = reshape(Rascot,[],1);
Z2Do = reshape(Zascot,[],1);
[inside,outside] = inpolygon(R2Do,Z2Do,RSURF(:,end), ZSURF(:,end));
in = inside & ! outside;
R2D = R2Do(in);
Z2D = Z2Do(in);
BMVOLo = (2*pi*( max(diff(Z2Do))* max(diff(R2Do)))*R2Do);
BMVOL = BMVOLo(in);
for i=1: length(P)-1
  A_D_NBI(i) = (P(i) + P(i+1))/2;
endfor
for i=1: length(E)-1
  E_D_NBI(i) = (E(i) + E(i+1))/2;
endfor

A_D_NBI = A_D_NBI';
E_D_NBI = E_D_NBI';
NE_D_NBI = length(E_D_NBI);
NA_D_NBI = length(A_D_NBI);
DEN = E_D_NBI(2) - E_D_NBI(1);
DPA = A_D_NBI(2) - A_D_NBI(1);
a=squeeze(sum(DEMO.distributions.rzPitchEdist.ordinate,7));
dist4adress = reshape(permute(a, [4,3,2,1]), [20, 40, 1250]);
F_D_NBIa = (dist4adress(:,:,in)./6.242e18)/5e5;
for i=1:length(BMVOL)
  F_D_NBI(:,:,i) = fliplr(F_D_NBIa(:,:,i));
end
NTOT_D_NBI = sum((BMVOL.*squeeze(sum(squeeze(sum(F_D_NBI,1)),1))')).*(DPA*DEN)/2
X2Do = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2Do, Z2Do);
X2D = X2Do(in);

TIMEfi = 5;
DT_AVG = 1;
nzones = length(ND);
N_2DZONES = length(BMVOL);
nsjdotb = -1;
nsnccw = -1;

FIDTGLOBAL = zeros(NE_D_NBI,NA_D_NBI);


for i=1:length(BMVOL)
  a = F_D_NBI(:,:,i).*BMVOL(i);
  FIDTGLOBAL = FIDTGLOBAL+a;
end 

[X Y] = meshgrid (linspace(min(Rhisto), max(Rhisto), length(Rhisto)),...
 linspace(min(Zhisto), max(Zhisto), length(Zhisto)));

w = sum(squeeze(sum(F_D_NBI,1)),1)*DPA*DEN/2;
FIDD = griddata(R2D, Z2D, w', X, Y, 'linear');
u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u


TOTRZ_TRANSP = 2*pi*sum(sum(FIDD.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))



FIDTP = sum(FIDTGLOBAL,1).*DEN;
FIDTE = sum(FIDTGLOBAL,2).*DPA;

figure(1)
subplot(1,3,1)
pcolor(-A_D_NBI, E_D_NBI/1000, FIDTGLOBAL)
title('FI density ev^{-1}', 'fontsize',14)
shading interp
xlabel('pitch ','fontsize', 14)
ylabel('energy (keV)','fontsize', 14)
shading interp
xlim([-1 1])
%ylim([0 70000])
hold on;
colormap(RBW);

colorbar('fontsize', 14)
set (gca, 'fontsize', 14)


subplot(1,3,2)
plot(-A_D_NBI, FIDTP,'r','linewidth',2)
xlabel('pitch ','fontsize', 14)
ylabel('dFID(P)/dP density','fontsize', 14)


subplot(1,3,3)
plot(E_D_NBI/1000, FIDTE,'r','linewidth',2)
xlabel('energy (keV)','fontsize', 14)
ylabel('dFID(E)/dE density','fontsize', 14)


figure(2)
pcolor(X/100,Y/100, FIDD)
hold on
plot(RSURF(:,end)./100, ZSURF(:,end)./100, 'linewidth',2,'k')
hold on
title('FI density m^{-3}', 'fontsize',14)
xlabel('major radius R (m)','fontsize', 14)
ylabel('Z (m)','fontsize', 14)
shading flat
axis equal
hold on;
colormap(RBW);
colorbar('fontsize', 14)
%caxis([0 4e12]);
set (gca, 'fontsize', 14)



figure(3)
subplot(1,3,1)
pcolor(R_grid, Z_grid, BRRZ')
hold on
plot(RSURF(:,end)./100, ZSURF(:,end)./100,'k')
title('B_R (T)')
shading flat
xlim([min(R_grid), max(R_grid)])
ylim([min(Z_grid), max(Z_grid)])
colorbar
xlabel('major radius (m)')
ylabel('Z (m)')


subplot(1,3,2)
pcolor(R_grid, Z_grid, BZRZ')
hold on
plot(RSURF(:,end)./100, ZSURF(:,end)./100,'k')
title('B_Z (T)')
shading flat
xlim([min(R_grid), max(R_grid)])
ylim([min(Z_grid), max(Z_grid)])
colorbar
xlabel('major radius (m)')
ylabel('Z (m)')


subplot(1,3,3)
pcolor(R_grid, Z_grid, BphiRZ')
hold on
plot(RSURF(:,end)./100, ZSURF(:,end)./100,'k')
title('B_{\phi} (T)')
shading flat
xlim([min(R_grid), max(R_grid)])
ylim([min(Z_grid), max(Z_grid)])
colorbar
xlabel('major radius (m)')
ylabel('Z (m)')

close all;

#I am adding these because they need to run dT properly


nccreate('DEMO_fi.cdf','NE_T_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','NA_T_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','NE_T_FUSN', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','NA_T_FUSN', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','E_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI)});
nccreate('DEMO_fi.cdf','A_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00040',length(A_D_NBI)});
nccreate('DEMO_fi.cdf','E_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI)});
nccreate('DEMO_fi.cdf','A_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00040',length(A_D_NBI)});

nccreate('DEMO_fi.cdf','F_T_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI),'dim_00040',length(A_D_NBI),'dim_00840',length(BMVOL)});
nccreate('DEMO_fi.cdf','F_T_FUSN', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI),'dim_00040',length(A_D_NBI),'dim_00840',length(BMVOL)});
nccreate('DEMO_fi.cdf','NTOT_T_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('DEMO_fi.cdf','NTOT_T_FUSN', 'Format', '64bit', 'Datatype', 'double');

nccreate('DEMO_fi.cdf','NE_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','NA_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','E_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI)});
nccreate('DEMO_fi.cdf','A_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00040',length(A_D_NBI)});
nccreate('DEMO_fi.cdf','F_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00020',length(E_D_NBI),'dim_00040',length(A_D_NBI),'dim_01250',length(BMVOL)});
nccreate('DEMO_fi.cdf','NTOT_D_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('DEMO_fi.cdf','NTHSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','R2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00699', length(BMVOL)});
nccreate('DEMO_fi.cdf','Z2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00699', length(BMVOL)});
nccreate('DEMO_fi.cdf','X2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00699', length(BMVOL)});
nccreate('DEMO_fi.cdf','BMVOL', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00699', length(BMVOL)});
nccreate('DEMO_fi.cdf','RSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00526',length(RSURF(:,1)),'dim_00050',length(RSURF(1,:))});
nccreate('DEMO_fi.cdf','ZSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00526',length(ZSURF(:,1)),'dim_00050',length(ZSURF(1,:))});
nccreate('DEMO_fi.cdf','TIME', 'Format', '64bit', 'Datatype', 'double');
nccreate('DEMO_fi.cdf','DT_AVG', 'Format', '64bit', 'Datatype', 'double');
nccreate('DEMO_fi.cdf','nzones', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','N_2DZONES', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','NXSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','nsjdotb', 'Format', '64bit', 'Datatype', 'int32');
nccreate('DEMO_fi.cdf','nsnccw', 'Format', '64bit', 'Datatype', 'int32');



ncwrite('DEMO_fi.cdf','NE_T_NBI',  NE_D_NBI);
ncwrite('DEMO_fi.cdf','NE_T_FUSN',  NE_D_NBI);

ncwrite('DEMO_fi.cdf','NA_T_NBI', NA_D_NBI);
ncwrite('DEMO_fi.cdf','NA_T_FUSN', NA_D_NBI);

ncwrite('DEMO_fi.cdf','E_T_NBI', E_D_NBI');
ncwrite('DEMO_fi.cdf','E_T_FUSN', E_D_NBI');

ncwrite('DEMO_fi.cdf','A_T_NBI', A_D_NBI');
ncwrite('DEMO_fi.cdf','A_T_FUSN', A_D_NBI');

ncwrite('DEMO_fi.cdf','F_T_NBI', F_D_NBI.*1e-20);
ncwrite('DEMO_fi.cdf','F_T_FUSN', F_D_NBI.*1e-20);

ncwrite('DEMO_fi.cdf','NTOT_T_NBI', 3.4705);
ncwrite('DEMO_fi.cdf','NTOT_T_FUSN', 3.4705);


ncwrite('DEMO_fi.cdf','NE_D_NBI', NE_D_NBI);
ncwrite('DEMO_fi.cdf','NA_D_NBI', NA_D_NBI);
ncwrite('DEMO_fi.cdf','E_D_NBI', E_D_NBI');
ncwrite('DEMO_fi.cdf','A_D_NBI', A_D_NBI');
ncwrite('DEMO_fi.cdf','NTHSURF', NTHSURF);
ncwrite('DEMO_fi.cdf','R2D', R2D);
ncwrite('DEMO_fi.cdf','Z2D', Z2D);
ncwrite('DEMO_fi.cdf','X2D', X2D);
ncwrite('DEMO_fi.cdf','BMVOL', BMVOL);
ncwrite('DEMO_fi.cdf','RSURF', RSURF);
ncwrite('DEMO_fi.cdf','ZSURF', ZSURF);
ncwrite('DEMO_fi.cdf','F_D_NBI', F_D_NBI);
ncwrite('DEMO_fi.cdf','NTOT_D_NBI', NTOT_D_NBI);
ncwrite('DEMO_fi.cdf','TIME', TIMEfi);
ncwrite('DEMO_fi.cdf','DT_AVG',  DT_AVG);
ncwrite('DEMO_fi.cdf','nzones', nzones);
ncwrite('DEMO_fi.cdf','N_2DZONES', N_2DZONES);
ncwrite('DEMO_fi.cdf','NXSURF', NXSURF);
ncwrite('DEMO_fi.cdf','nsjdotb', nsjdotb);
ncwrite('DEMO_fi.cdf','nsnccw', nsnccw);


return

nccreate('DEMO.CDF','TIME', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME', length(TIME)});
nccreate('DEMO.CDF','TIME3', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME3', length(TIME3)});
nccreate('DEMO.CDF','X', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', length(X),'TIME3', length(TIME3)});
nccreate('DEMO.CDF','OMEGA', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', length(X),'TIME3', length(TIME3)});
nccreate('DEMO.CDF','ND', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', length(X),'TIME3', length(TIME3)});
nccreate('DEMO.CDF','NT', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', length(X),'TIME3', length(TIME3)});
nccreate('DEMO.CDF','TI', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', length(X),'TIME3', length(TIME3)});

ncwrite('DEMO.CDF','TIME', TIME);
ncwrite('DEMO.CDF','TIME3', TIME3);
ncwrite('DEMO.CDF','X', X);
ncwrite('DEMO.CDF','OMEGA', OMEGA);
ncwrite('DEMO.CDF','ND', ND);
ncwrite('DEMO.CDF','NT', NT);
ncwrite('DEMO.CDF','TI', TI);

nccreate('DEMO_ps_state.cdf','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', length(R_grid)});
nccreate('DEMO_ps_state.cdf','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', length(Z_grid)});
nccreate('DEMO_ps_state.cdf','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(R_grid),'dim_nz', length(Z_grid)});
nccreate('DEMO_ps_state.cdf','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(R_grid),'dim_nz', length(Z_grid)});
nccreate('DEMO_ps_state.cdf','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(R_grid),'dim_nz', length(Z_grid)});


ncwrite('DEMO_ps_state.cdf','R_grid', R_grid);
ncwrite('DEMO_ps_state.cdf','Z_grid', Z_grid);
ncwrite('DEMO_ps_state.cdf','BRRZ', BRRZ);
ncwrite('DEMO_ps_state.cdf','BphiRZ', BphiRZ);
ncwrite('DEMO_ps_state.cdf','BZRZ', BZRZ);