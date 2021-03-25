%%%%%%%%%%%%%%%%%5 HERE I READ SOME VARIABLE FROM TRANSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fast = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29210/O22/29210O22_fi_6.cdf';
%state = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_ps_ts1_state.cdf';
big = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29210/O22/29210O22.CDF';
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

%{
%PLASMA STATE
R_grid = ncread(state, 'R_grid');
Z_grid = ncread(state, 'Z_grid');
BRRZ = ncread(state, 'BRRZ');
BphiRZ = ncread(state, 'BphiRZ');
BZRZ = ncread(state, 'BZRZ');
%}
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


RSURF = ncread(fast, 'RSURF');
ZSURF = ncread(fast, 'ZSURF');
F_D_NBI = ncread(fast, 'F_D_NBI');

NTOT_D_NBI = ncread(fast, 'NTOT_D_NBI');

TIMEfi = ncread(fast, 'TIME');
DT_AVG = ncread(fast, 'DT_AVG');
nzones = ncread(fast, 'nzones');
N_2DZONES = ncread(fast,'N_2DZONES');
NXSURF = ncread(fast,'NXSURF');
#nsjdotb = ncread(fast,'nsjdotb');
#nsnccw = ncread(fast,'nsnccw');

%%%%%%%%%%%%%%%%%5 POLOIDAL FLUX COORDINATES FOR ASCOT5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%5

XB = ncread(big, 'XB');
PLFLX = ncread(big, 'PLFLX');
rho_pol = [0;PLFLX(:,56)];
rho_pol = sqrt(rho_pol/rho_pol(end));
X_A = interp1([0;XB(:,56)], rho_pol, Xrho(:,56));
X_A = X_A.*ones(60,58);

newdens = interp1(X_A(:,56),ND(:,56), Xrho(:,56));
newtemp = interp1(X_A(:,56),TI(:,56), Xrho(:,56));
newomega = interp1(X_A(:,56),OMEGA(:,56), Xrho(:,56));
electrontemp = interp1(X_A(:,56),TE(:,56), Xrho(:,56));
electrondens = interp1(X_A(:,56),NE(:,56), Xrho(:,56));


newdens(1) = ND(1,56);
newtemp(1) = TI(1,56);
newomega(1) = OMEGA(1,56);
ND_A  = newdens.*ones(60,58);
TI_A  = newtemp.*ones(60,58);
OMEGA_A  = newomega.*ones(60,58);



%%%%%%%%%%%%%%%%%5 HERE I TAKE INTO ACCOUNT ONLY POINTS INSIDE THE LCFS AND I CALCULATE SOME OTHER VARIABLES FOR DRESS %%%%%%%%

addpath('/home/andrea/Documents/MASTOrbit');
#[EQD] = read_eqdsk('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_matlab/g29210.00255',0);

R = linspace(0.21, 1.43, 62);
Z = linspace(-1.01, 1.21, 112);


for i=1: length(R)-1
  Rhisto(i) = (R(i) + R(i+1))/2;
end
Rhisto = Rhisto.*100;%

for i=1: length(Z)-1
  Zhisto(i) = (Z(i) + Z(i+1))/2;
end
Zhisto = Zhisto.*100;


% Interpolates the FID on the grid for plotting
DEN = E_D_NBI(2) - E_D_NBI(1);
DPA = A_D_NBI(2) - A_D_NBI(1);


[Rascot, Zascot] = meshgrid(Rhisto, Zhisto);
R2D_Ao = reshape(Rascot,[],1);
Z2D_Ao = reshape(Zascot,[],1);
[inside,outside] = inpolygon(R2D_Ao,Z2D_Ao,RSURF(:,end), ZSURF(:,end));
in = inside & ! outside;
R2D_A = R2D_Ao(in);
Z2D_A = Z2D_Ao(in);
X2D_Ao = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2D_Ao, Z2D_Ao);
X2D_A = X2D_Ao(in);
#X2D_A = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2D_A, Z2D_A);
BMVOL_Ao = (2*pi*( max(diff(Z2D_Ao))* max(diff(R2D_Ao)))*R2D_Ao);
BMVOL_A = BMVOL_Ao(in);



%%%%%%%%%%%%% ASCOT 4 %%%%%%%%%%%%%%%%%%%%%
ascot4_SS = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29210/SS/SS_D.h5');
dist4a_SS = squeeze(ascot4_SS.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adresso_SS = reshape(permute(dist4a_SS, [4,3,2,1]), [75, 50, 6771]);
ascot4_SW = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29210/SW/SW_D.h5');
dist4a_SW = squeeze(ascot4_SW.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adresso_SW = reshape(permute(dist4a_SW, [4,3,2,1]), [75, 50, 6771]);
dist4adresso = dist4adresso_SW + dist4adresso_SS;


dist4adressa = dist4adresso(:,:,in)/569173;
dist4adress = dist4adressa;
for i=1:4362
  dist4adress(:,:,i) = fliplr(dist4adressa(:,:,i));
end


clear dist4adresso;

FIDTGLOBAL = zeros(75,50);
FIDA4GLOBAL = zeros(75,50);

for i=1:length(BMVOL)
  a = F_D_NBI(:,:,i).*BMVOL(i);
  FIDTGLOBAL = FIDTGLOBAL+a;
end 


for i=1:length(BMVOL_A)
  t = dist4adress(:,:,i).*BMVOL_A(i);
  FIDA4GLOBAL = FIDA4GLOBAL+t;
end 

pitch4 = ascot4_SW.distributions.rzPitchEdist.abscissae.dim3;
energy4 = ascot4_SW.distributions.rzPitchEdist.abscissae.dim4.*6.242e18;
for i=1: length(pitch4)-1
  p4(i) = (pitch4(i) + pitch4(i+1))/2;
endfor

for i=1: length(energy4)-1
  e4(i) = (energy4(i) + energy4(i+1))/2;
endfor

[X Y] = meshgrid (linspace(min(Rhisto), max(Rhisto), length(Rhisto)),...
 linspace(min(Zhisto), max(Zhisto), length(Zhisto)));

w = sum(squeeze(sum(F_D_NBI,1)),1)*DPA*DEN/2;
h = sum(squeeze(sum(dist4adress,1)),1)*(p4(2)-p4(1)).*(e4(2)-e4(1))/2;
FIDD = griddata(R2D, Z2D, w', X, Y, 'linear');
FIDDA4 = griddata(R2D_A, Z2D_A, h', X, Y, 'linear');

u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u
u = find(isnan(FIDDA4) == 1);
FIDDA4(u) = 0;
clear u




figure(1)
subplot(1,3,1)
pcolor(X,Y, FIDD)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
#plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('TRANSP , FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0 3e12]);
set (gca, 'fontsize', 14)

subplot(1,3,2)
pcolor(X,Y, FIDDA4)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
#plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('ASCOT 4 in FO, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0 3e12]);
set (gca, 'fontsize', 14)

subplot(1,3,3)
pcolor(X,Y, 100.*(1-FIDD./FIDDA4))
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
#plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('Relative difference', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([-100 100]);
set (gca, 'fontsize', 14)




FIDTP = sum(FIDTGLOBAL,1).*DEN;
FIDTE = sum(FIDTGLOBAL,2).*DPA;
FIDA4P = sum(FIDA4GLOBAL,1).*(e4(2)-e4(1));
FIDA4E = sum(FIDA4GLOBAL,2).*(p4(2)-p4(1));

figure(100)
subplot(1,3,1)
pcolor(-A_D_NBI, E_D_NBI,  FIDTGLOBAL)
shading interp
title('TRANSP w/o FLR, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 75000])
hold on;
colormap(RBW);
caxis([0 4e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)

subplot(1,3,2)
pcolor(-p4, e4, FIDA4GLOBAL)
shading interp
title('ASCOT4 in GC, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 75000])
hold on;
colormap(RBW);
caxis([0 4e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)

subplot(1,3,3)
pcolor(-A_D_NBI, E_D_NBI, 100.*(1-FIDTGLOBAL./FIDA4GLOBAL))
shading interp
title('Relative difference', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 75000])
hold on;
colormap(RBW);
caxis([-100 100]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)




figure(4)
subplot(2,2,1)
plot(A_D_NBI, FIDTP,'r','linewidth',2,...
p4, FIDA4P,'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('pitch')
ylabel('dFID(P)/dP density')


subplot(2,2,3)
plot(E_D_NBI, FIDTE,'r','linewidth',2,...
e4, FIDA4E,'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')

subplot(2,2,3)
plot(A_D_NBI, FIDTP./sum(FIDTP),'r','linewidth',2,...
p4, FIDA4P./sum(FIDA4P),'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('pitch')
ylabel('dFID(P)/dP density')

subplot(2,2,4)
plot(E_D_NBI, FIDTE./sum(FIDTE),'r','linewidth',2,...
e4, FIDA4E./sum(FIDA4E),'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')





TOTEP_TRANSP= sum(FIDTGLOBAL(:)).*DPA.*DEN/2
TOTEP_ASCOT4 = sum(FIDA4GLOBAL(:)).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2

TOTRZ_TRANSP = 2*pi*sum(sum(FIDD.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))
TOTRZ_ASCOT4 = 2*pi*sum(sum(FIDDA4.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))

ASCOT4 = sum((BMVOL_A.*squeeze(sum(squeeze(sum(dist4adress,1)),1))')).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2
TRANSP = sum((BMVOL  .*squeeze(sum(squeeze(sum(F_D_NBI,1)),1))')).*(DPA*DEN)/2


%%%%%%%%%%%%%%%% FAST IONS FILE FROM ASCOT %%%%%%%%%%%%%%%
return

nccreate('29210_fi_1.cdf','NE_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','NA_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','E_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4)});
nccreate('29210_fi_1.cdf','A_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(p4)});
nccreate('29210_fi_1.cdf','NTHSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','R2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04362', length(BMVOL_A)});
nccreate('29210_fi_1.cdf','Z2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04362', length(BMVOL_A)});
nccreate('29210_fi_1.cdf','X2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04362', length(BMVOL_A)});
nccreate('29210_fi_1.cdf','BMVOL', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04362', length(BMVOL_A)});
nccreate('29210_fi_1.cdf','RSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_0041',41});
nccreate('29210_fi_1.cdf','ZSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_0041',41});
nccreate('29210_fi_1.cdf','F_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4),'dim_00050',length(p4),'dim_04362',length(BMVOL_A)});
nccreate('29210_fi_1.cdf','NTOT_D_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('29210_fi_1.cdf','TIME', 'Format', '64bit', 'Datatype', 'double');
nccreate('29210_fi_1.cdf','DT_AVG', 'Format', '64bit', 'Datatype', 'double');
nccreate('29210_fi_1.cdf','nzones', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','N_2DZONES', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','NXSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','nsjdotb', 'Format', '64bit', 'Datatype', 'int32');
nccreate('29210_fi_1.cdf','nsnccw', 'Format', '64bit', 'Datatype', 'int32');

ncwrite('29210_fi_1.cdf','NE_D_NBI', length(e4));
ncwrite('29210_fi_1.cdf','NA_D_NBI', length(p4));
ncwrite('29210_fi_1.cdf','E_D_NBI', e4');
ncwrite('29210_fi_1.cdf','A_D_NBI', p4');
ncwrite('29210_fi_1.cdf','NTHSURF', NTHSURF);
ncwrite('29210_fi_1.cdf','R2D', R2D_A);
ncwrite('29210_fi_1.cdf','Z2D', Z2D_A);
ncwrite('29210_fi_1.cdf','X2D', X2D_A);
ncwrite('29210_fi_1.cdf','BMVOL', BMVOL_A);
ncwrite('29210_fi_1.cdf','RSURF', RSURF);
ncwrite('29210_fi_1.cdf','ZSURF', ZSURF);
ncwrite('29210_fi_1.cdf','F_D_NBI', dist4adress);
ncwrite('29210_fi_1.cdf','NTOT_D_NBI', ASCOT4);
ncwrite('29210_fi_1.cdf','TIME', TIMEfi);
ncwrite('29210_fi_1.cdf','DT_AVG',  DT_AVG);
ncwrite('29210_fi_1.cdf','nzones', 60);
ncwrite('29210_fi_1.cdf','N_2DZONES', length(BMVOL_A));
ncwrite('29210_fi_1.cdf','NXSURF', NXSURF);
ncwrite('29210_fi_1.cdf','nsjdotb', -1);
ncwrite('29210_fi_1.cdf','nsnccw', -1);


%%%%%%%%%%%%%%%% CDF FILE %%%%%%%%%%%%%%%

nccreate('29210.CDF','TIME', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME', length(TIME)});
nccreate('29210.CDF','TIME3', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME3', length(TIME3)});
nccreate('29210.CDF','X', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 163,'TIME3', length(TIME3)});
nccreate('29210.CDF','OMEGA', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 163,'TIME3', length(TIME3)});
nccreate('29210.CDF','ND', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 163,'TIME3', length(TIME3)});
nccreate('29210.CDF','TI', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X',163,'TIME3', length(TIME3)});

ncwrite('29210.CDF','TIME', TIME);
ncwrite('29210.CDF','TIME3', TIME3);
ncwrite('29210.CDF','X', Xrho);
ncwrite('29210.CDF','OMEGA', OMEGA_A);
ncwrite('29210.CDF','ND', ND_A);
ncwrite('29210.CDF','TI', TI_A);

%{
%%%%%%%%%%%%%%%% PLASMA STATE %%%%%%%%%%%%%%%
R_grid = ncread(state, 'R_grid');
Z_grid = ncread(state, 'Z_grid');
BRRZ = ncread(state, 'BRRZ');
BphiRZ = ncread(state, 'BphiRZ');
BZRZ = ncread(state, 'BZRZ');

nccreate('29210_ps_ts1_state.cdf','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', 105});
nccreate('29210_ps_ts1_state.cdf','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', 165});
nccreate('29210_ps_ts1_state.cdf','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('29210_ps_ts1_state.cdf','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('29210_ps_ts1_state.cdf','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});

ncwrite('29210_ps_ts1_state.cdf','R_grid', R_grid);
ncwrite('29210_ps_ts1_state.cdf','Z_grid', Z_grid);
ncwrite('29210_ps_ts1_state.cdf','BRRZ', BRRZ);
ncwrite('29210_ps_ts1_state.cdf','BphiRZ', BphiRZ);
ncwrite('29210_ps_ts1_state.cdf','BZRZ', BZRZ);
%}


nccreate('29210_ps_ts1_state.cdf','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R)});
nccreate('29210_ps_ts1_state.cdf','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', length(EQD.Z)});
nccreate('29210_ps_ts1_state.cdf','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});
nccreate('29210_ps_ts1_state.cdf','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});
nccreate('29210_ps_ts1_state.cdf','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});

ncwrite('29210_ps_ts1_state.cdf','R_grid', EQD.R');
ncwrite('29210_ps_ts1_state.cdf','Z_grid', EQD.Z');
ncwrite('29210_ps_ts1_state.cdf','BRRZ', EQD.B.POL.R);
ncwrite('29210_ps_ts1_state.cdf','BphiRZ', EQD.B.TOR);
ncwrite('29210_ps_ts1_state.cdf','BZRZ', EQD.B.POL.Z);


  

