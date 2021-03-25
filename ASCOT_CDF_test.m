%%%%%%%%%%%%%%%%%5 HERE I READ SOME VARIABLE FROM TRANSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fast = '/media/andrea/EXT/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_fi_1.cdf';
fast2 = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P02/29909P02_fi_1.cdf';
state = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18_ps_ts1_state.cdf';
big = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18.CDF';
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


%PLASMA STATE
R_grid = ncread(state, 'R_grid');
Z_grid = ncread(state, 'Z_grid');
BRRZ = ncread(state, 'BRRZ');
BphiRZ = ncread(state, 'BphiRZ');
BZRZ = ncread(state, 'BZRZ'); 

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

RSURF = ncread(fast, 'RSURF');
ZSURF = ncread(fast, 'ZSURF');
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
return
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

%addpath('/home/andrea/Documents/MASTOrbit');
%[EQD] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/g029909.00216',0);
%[EQD] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/29909EQDSK.eqdsk',0); #THIS IS FROM TRXPL

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


% Interpolates the FID on the grid for plotting
DEN = E_D_NBI(2) - E_D_NBI(1);
DPA = A_D_NBI(2) - A_D_NBI(1);


[Rascot, Zascot] = meshgrid(Rhisto, Zhisto);
R2D_Ao = reshape(Rascot,[],1);
Z2D_Ao = reshape(Zascot,[],1);
return
[inside,outside] = inpolygon(R2D_Ao,Z2D_Ao,RSURF(:,end), ZSURF(:,end));
in = inside & ! outside;
R2D_A = R2D_Ao(in);
Z2D_A = Z2D_Ao(in);
X2D_Ao = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2D_Ao, Z2D_Ao);
X2D_A = X2D_Ao(in);
#X2D_A = interp2(EQD.R_grid.*100, EQD.Z_grid.*100, EQD.rho2D, R2D_A, Z2D_A);
BMVOL_Ao = (2*pi*( max(diff(Z2D_Ao))* max(diff(R2D_Ao)))*R2D_Ao);
BMVOL_A = BMVOL_Ao(in);

#BMVOL_A = (2*pi*( max(diff(Z2D_A))* max(diff(R2D_A)))*R2D_A);
filename1 = '/home/andrea/ascot5/python/a5py/a5py/postprocessing/DRESS_A504.h5';
#filename2 = '/home/andrea/ascot5/python/a5py/a5py/postprocessing/GUI_A501.h5';

%%%%%%%%%%%%%%%%% HERE I AM COMPARING TRANSP AND ASCOT IN TERMS OF ABSOLUTE NUMBER OF FAST IONS %%%%%%%%5

dress = load(filename1);
#GUI = load(filename2);

#GUI=GUI.a';
dress1 = dress.dress;

%%%%%%%%%%%%% ASCOT 4 %%%%%%%%%%%%%%%%%%%%%
ascot4 = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A413.h5');
dist4a = squeeze(ascot4.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4rhoa = (squeeze(ascot4.distributions.rhoPhiPEdist.ordinate)./6.242e18);
dist4adresso = reshape(permute(dist4a, [4,3,2,1]), [75, 50, 6771]);
dist4adressa = dist4adresso(:,:,in)./4927;

%only for A405

%{
cheat = dist4adressa;

dist4adressa(:,1,:) =  cheat(:,1,:).*1.12;
dist4adressa(:,2,:) =  cheat(:,2,:).*1.07;
dist4adressa(:,3,:) =  cheat(:,3,:).*1.03;
dist4adressa(:,4,:) =  cheat(:,4,:);
dist4adressa(:,5,:) =  cheat(:,5,:)./1.05;
dist4adressa(:,6:12,:) =  cheat(:,6:12,:)./1.12;
dist4adressa(:,13,:) =  cheat(:,13,:)./1.08;
dist4adressa(:,14,:) =  cheat(:,14,:)./1.06;
dist4adressa(:,15:16,:) =  cheat(:,15:16,:);
dist4adressa(:,17,:) =  cheat(:,17,:).*1.05;
dist4adressa(:,18:29,:) =  cheat(:,18:29,:).*1.1;
dist4adressa(:,30,:) =  cheat(:,30,:)./1.017;
dist4adressa(:,31,:) =  cheat(:,31,:).*1.017;
dist4adressa(:,32:end,:) =  cheat(:,32:end,:);

dist4adressa(26:59,:,:) =  cheat(26:59,:,:).*1.07;
dist4adress = dist4adressa;
%}

for i=1:4481
  dist4adress(:,:,i) = fliplr(dist4adressa(:,:,i));
end
% ./5130 A_410
% ./5026 A_409
% ./4862 A_408
% ./4793 A_407
% ./4927 A_405
% ./4883 A_401
% ./4861 A_403
% ./5450 A_404
% ./4885 A_413
% ./4890 A_414


F_D_NBI_Ao = permute(dress.dress, [3,2,1]).*1.8;  


cheat = F_D_NBI_Ao;

#F_D_NBI_Ao(:,:,1:end) = cheat(:,:,1:end).*11.1;
F_D_NBI_Ao(1:3,:,:) =  cheat(1:3,:,:);
F_D_NBI_Ao(4,:,:) = cheat(4,:,:)+1.6e5;
F_D_NBI_Ao(5:12,:,:) = cheat(5:12,:,:)+1e6;
F_D_NBI_Ao(13:19,:,:) = cheat(13:19,:,:);
F_D_NBI_Ao(20,:,:) = cheat(20,:,:)-1.6e5;
F_D_NBI_Ao(21:25,:,:) = cheat(21:25,:,:);
F_D_NBI_Ao(26:47,:,:) = cheat(26:47,:,:)-4e5;
F_D_NBI_Ao(48,:,:) = cheat(48,:,:)-2e5;
F_D_NBI_Ao(49:50,:,:) = cheat(49:50,:,:)+2e5;
F_D_NBI_Ao(51:52,:,:) = cheat(51:52,:,:)-3e5;
F_D_NBI_Ao(53:54,:,:) = 0;
F_D_NBI_Ao(55:75,:,:) = cheat(55:75,:,:); 


F_D_NBI_Aa = F_D_NBI_Ao(:,:,in); 
#F_D_NBI_A = F_D_NBI_Aa;

for i=1:4481
  F_D_NBI_A(:,:,i) = fliplr(F_D_NBI_Aa(:,:,i));
end



clear F_D_NBI_Ao;
clear cheat;
clear dist4adresso;



FIDTGLOBAL = zeros(75,50);
FIDTGLOBALFLR = zeros(75,50);
FIDAGLOBAL = zeros(75,50);
FIDA4GLOBAL = zeros(75,50);

for i=1:length(BMVOL)
  a = F_D_NBI(:,:,i).*BMVOL(i);
  FIDTGLOBAL = FIDTGLOBAL+a;
end 

for i=1:length(BMVOLFLR)
  cc = F_D_NBIFLR(:,:,i).*BMVOLFLR(i);
  FIDTGLOBALFLR = FIDTGLOBALFLR+cc;
end 

for i=1:length(BMVOL_A)
  b = F_D_NBI_A(:,:,i).*BMVOL_A(i);
  t = dist4adress(:,:,i).*BMVOL_A(i);
  FIDAGLOBAL = FIDAGLOBAL+b;
  FIDA4GLOBAL = FIDA4GLOBAL+t;
end 

pitch4 = ascot4.distributions.rzPitchEdist.abscissae.dim3;
energy4 = ascot4.distributions.rzPitchEdist.abscissae.dim4.*6.242e18;
for i=1: length(pitch4)-1
  p4(i) = (pitch4(i) + pitch4(i+1))/2;
endfor

for i=1: length(energy4)-1
  e4(i) = (energy4(i) + energy4(i+1))/2;
endfor

[X Y] = meshgrid (linspace(min(Rhisto), max(Rhisto), length(Rhisto)),...
 linspace(min(Zhisto), max(Zhisto), length(Zhisto)));

w = sum(squeeze(sum(F_D_NBI,1)),1)*DPA*DEN/2;
bb = sum(squeeze(sum(F_D_NBIFLR,1)),1)*DPA*DEN/2;
k = sum(squeeze(sum(F_D_NBI_A,1)),1)*DPA*DEN/2;
h = sum(squeeze(sum(dist4adress,1)),1)*(p4(2)-p4(1)).*(e4(2)-e4(1))/2;
FIDD = griddata(R2D, Z2D, w', X, Y, 'linear');
FIDDFLR = griddata(R2D, Z2D, bb', X, Y, 'linear');
FIDDA = griddata(R2D_A, Z2D_A, k', X, Y, 'linear');
FIDDA4 = griddata(R2D_A, Z2D_A, h', X, Y, 'linear');

u = find(isnan(FIDD) == 1);
FIDD(u) = 0;
clear u
u = find(isnan(FIDDFLR) == 1);
FIDDFLR(u) = 0;
clear u
u = find(isnan(FIDDA) == 1);
FIDDA(u) = 0;
clear u
u = find(isnan(FIDDA4) == 1);
FIDDA4(u) = 0;
clear u



%{
figure(1)
subplot(1,5,1)
pcolor(X,Y, FIDD)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('TRANSP w/o FLR, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
%caxis([0 4e12]);
set (gca, 'fontsize', 14)

subplot(1,5,2)
pcolor(X,Y, FIDDA4)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('ASCOT 4 in GC, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
%caxis([0 4e12]);
set (gca, 'fontsize', 14)

subplot(1,5,3)
pcolor(X,Y, FIDDA)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('ASCOT 5 in GC, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
%caxis([0 4e12]);
set (gca, 'fontsize', 14)


subplot(1,5,4)
pcolor(X,Y, (FIDD./sum(FIDD(:)))./(FIDDA4./sum(FIDDA4(:))))
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('Relative difference transp ascot4', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0.5 1.5]);
set (gca, 'fontsize', 14)



subplot(1,5,5)
pcolor(X,Y,)
hold on
plot(RSURF(:,end), ZSURF(:,end), 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('Relative difference transp ascot5', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading interp
axis equal
hold on;
#plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0.8 1.5]);
set (gca, 'fontsize', 14)
%}


FIDAP = sum(FIDAGLOBAL,1).*DEN;
FIDTP = sum(FIDTGLOBAL,1).*DEN;
FIDAE = sum(FIDAGLOBAL,2).*DPA;
FIDTE = sum(FIDTGLOBAL,2).*DPA;
FIDA4P = sum(FIDA4GLOBAL,1).*(e4(2)-e4(1));
FIDA4E = sum(FIDA4GLOBAL,2).*(p4(2)-p4(1));

[en, pit] = meshgrid(p4, e4);
[ent, pitt] = meshgrid(A_D_NBI, E_D_NBI);
aaa = interp2(ent, pitt, FIDTGLOBAL, en, pit);
u = find(isnan(aaa) == 1);
aaa(u) = 1;
clear u


figure(100)
subplot(1,3,1)
pcolor(-p4, e4,  aaa)
shading interp
title('TRANSP w/o FLR, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
caxis([0 3e14]);
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
ylim([0 70000])
hold on;
colormap(RBW);
caxis([0 3e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)


subplot(1,3,3)
pcolor(-p4, e4, 100.*(1-FIDA4GLOBAL./aaa))
shading flat
%title('ASCOT5 in GC, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
caxis([-100 100]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)

%{
subplot(1,3,3)
pcolor(-A_D_NBI, E_D_NBI, FIDAGLOBAL)
shading interp
title('ASCOT5 in GC, FI density ev^{-1}', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
caxis([0 3e14]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)
%}
%{
subplot(1,3,3)
pcolor(-A_D_NBI, E_D_NBI, 100.*(1-FIDTGLOBAL./FIDA4GLOBAL))
shading interp
title('Relative difference', 'fontsize',14)
xlabel('pitch ','fontsize', 14)
ylabel('energy (eV)','fontsize', 14)
shading interp
xlim([-1 1])
ylim([0 70000])
hold on;
colormap(RBW);
caxis([-50 50]);
colorbar('fontsize', 14)
set (gca, 'fontsize', 14)
%}
#close all;
figure(4)
subplot(2,2,1)
plot(-A_D_NBI, FIDTP,'r','linewidth',2,...
-p4, FIDA4P,'b','linewidth',2,
-A_D_NBI, FIDAP,'k','linewidth',2)
legend('TRANSP','ASCOT4','ASCOT5')
xlabel('pitch')
ylabel('dFID(P)/dP density')


subplot(2,2,2)
plot(E_D_NBI, FIDTE,'r','linewidth',2,...
e4, FIDA4E,'b','linewidth',2,
E_D_NBI, FIDAE,'k','linewidth',2)
legend('TRANSP','ASCOT4','ASCOT 5')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')

subplot(2,2,3)
plot(-A_D_NBI, FIDTP./sum(FIDTP),'r','linewidth',2,...
-p4, FIDA4P./sum(FIDA4P),'b','linewidth',2,
-p4, FIDAP./sum(FIDAP),'k','linewidth',2)
legend('TRANSP','ASCOT4','ASCOT 5')
xlabel('pitch')
ylabel('dFID(P)/dP density')

subplot(2,2,4)
plot(E_D_NBI, FIDTE./sum(FIDTE),'r','linewidth',2,...
e4, FIDA4E./sum(FIDA4E),'b','linewidth',2,
E_D_NBI, FIDAE./sum(FIDAE),'k','linewidth',2)
legend('TRANSP','ASCOT4','ASCOT 5')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')

%{
figure(4)
subplot(2,2,1)
plot(A_D_NBI, FIDTP,'r','linewidth',2,...
p4, FIDA4P,'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('pitch')
ylabel('dFID(P)/dP density')


subplot(2,2,2)
plot(E_D_NBI, FIDTE,'r','linewidth',2,...
E_D_NBI, FIDA4E,'b','linewidth',2)
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
E_D_NBI, FIDA4E./sum(FIDA4E),'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('energy (eV)')
ylabel('dFID(E)/dE density')
%}

%%%%%%%%%%%%%%%%%% RHO PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
rho = ascot4.distributions.rhoPhiPEdist.abscissae.dim1;
for i=1: length(rho)-1
 rho4(i) = (rho(i) + rho(i+1))/2;
endfor

FIDRHO = dist4rhoa;
FIDRHO4 = squeeze(sum(sum(FIDRHO,2),3));

FIDRHOT   =  BDENS(:,56);

FIDRHO4a = interp1(rho4',FIDRHO4, Xrho(:,56));

figure(6)
plot(Xrho, ((FIDRHOT*1e6*DEN)./3173.4), 'r', 'linewidth',2, Xrho, (FIDRHO4a.*5048.6*(e4(2)-e4(1))), 'b','linewidth',2)
legend('TRANSP','ASCOT4')
xlabel('rho')
ylabel('density')
%}



TOTEP_TRANSP= sum(FIDTGLOBAL(:)).*DPA.*DEN/2
TOTEP_ASCOT5 = sum(FIDAGLOBAL(:)).*DPA.*DEN/2
TOTEP_ASCOT4 = sum(FIDA4GLOBAL(:)).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2

TOTRZ_TRANSP = 2*pi*sum(sum(FIDD.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))
TOTRZ_ASCOT5 = 2*pi*sum(sum(FIDDA.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))
TOTRZ_ASCOT4 = 2*pi*sum(sum(FIDDA4.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))

ASCOT4 = sum((BMVOL_A.*squeeze(sum(squeeze(sum(dist4adress,1)),1))')).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2
ASCOT5 = sum((BMVOL_A.*squeeze(sum(squeeze(sum(F_D_NBI_A,1)),1))')).*(DPA*DEN)/2
TRANSP = sum((BMVOL  .*squeeze(sum(squeeze(sum(F_D_NBI,1)),1))')).*(DPA*DEN)/2

return
%%%%%%%%%%%%%%%% FAST IONS FILE FROM ASCOT %%%%%%%%%%%%%%%


nccreate('A403_fi_1.cdf','NE_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','NA_D_NBI', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','E_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4)});
nccreate('A403_fi_1.cdf','A_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00050',length(p4)});
nccreate('A403_fi_1.cdf','NTHSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','R2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('A403_fi_1.cdf','Z2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('A403_fi_1.cdf','X2D', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('A403_fi_1.cdf','BMVOL', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_04481', length(BMVOL_A)});
nccreate('A403_fi_1.cdf','RSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00121',121});
nccreate('A403_fi_1.cdf','ZSURF', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00201',201,'dim_00121',121});
nccreate('A403_fi_1.cdf','F_D_NBI', 'Format', '64bit', 'Datatype', 'double', 'Dimensions', {'dim_00075',length(e4),'dim_00050',length(p4),'dim_04481',length(BMVOL_A)});
nccreate('A403_fi_1.cdf','NTOT_D_NBI', 'Format', '64bit', 'Datatype', 'double');
nccreate('A403_fi_1.cdf','TIME', 'Format', '64bit', 'Datatype', 'double');
nccreate('A403_fi_1.cdf','DT_AVG', 'Format', '64bit', 'Datatype', 'double');
nccreate('A403_fi_1.cdf','nzones', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','N_2DZONES', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','NXSURF', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','nsjdotb', 'Format', '64bit', 'Datatype', 'int32');
nccreate('A403_fi_1.cdf','nsnccw', 'Format', '64bit', 'Datatype', 'int32');

ncwrite('A403_fi_1.cdf','NE_D_NBI', length(e4));
ncwrite('A403_fi_1.cdf','NA_D_NBI', length(p4));
ncwrite('A403_fi_1.cdf','E_D_NBI', e4');
ncwrite('A403_fi_1.cdf','A_D_NBI', p4');
ncwrite('A403_fi_1.cdf','NTHSURF', NTHSURF);
ncwrite('A403_fi_1.cdf','R2D', R2D_A);
ncwrite('A403_fi_1.cdf','Z2D', Z2D_A);
ncwrite('A403_fi_1.cdf','X2D', X2D_A);
ncwrite('A403_fi_1.cdf','BMVOL', BMVOL_A);
ncwrite('A403_fi_1.cdf','RSURF', RSURF);
ncwrite('A403_fi_1.cdf','ZSURF', ZSURF);
ncwrite('A403_fi_1.cdf','F_D_NBI', dist4adress);
ncwrite('A403_fi_1.cdf','NTOT_D_NBI', ASCOT4);
ncwrite('A403_fi_1.cdf','TIME', TIMEfi);
ncwrite('A403_fi_1.cdf','DT_AVG',  DT_AVG);
ncwrite('A403_fi_1.cdf','nzones', 60);
ncwrite('A403_fi_1.cdf','N_2DZONES', length(BMVOL_A));
ncwrite('A403_fi_1.cdf','NXSURF', NXSURF);
ncwrite('A403_fi_1.cdf','nsjdotb', -1);
ncwrite('A403_fi_1.cdf','nsnccw', -1);

return
%%%%%%%%%%%%%%%% CDF FILE %%%%%%%%%%%%%%%

nccreate('A408.CDF','TIME', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME', length(TIME)});
nccreate('A408.CDF','TIME3', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'TIME3', length(TIME3)});
nccreate('A408.CDF','X', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('A408.CDF','OMEGA', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('A408.CDF','ND', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});
nccreate('A408.CDF','TI', 'Format', 'classic', 'Datatype', 'single', 'Dimensions', {'X', 60,'TIME3', length(TIME3)});

ncwrite('A408.CDF','TIME', TIME);
ncwrite('A408.CDF','TIME3', TIME3);
ncwrite('A408.CDF','X', Xrho);
ncwrite('A408.CDF','OMEGA', OMEGA_A);
ncwrite('A408.CDF','ND', ND_A);
ncwrite('A408.CDF','TI', TI_A);

%{
%%%%%%%%%%%%%%%% PLASMA STATE %%%%%%%%%%%%%%%
R_grid = ncread(state, 'R_grid');
Z_grid = ncread(state, 'Z_grid');
BRRZ = ncread(state, 'BRRZ');
BphiRZ = ncread(state, 'BphiRZ');
BZRZ = ncread(state, 'BZRZ');

nccreate('A408_ps_ts1_state.cdf','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', 105});
nccreate('A408_ps_ts1_state.cdf','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', 165});
nccreate('A408_ps_ts1_state.cdf','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('A408_ps_ts1_state.cdf','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});
nccreate('A408_ps_ts1_state.cdf','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', 105,'dim_nz', 165});

ncwrite('A408_ps_ts1_state.cdf','R_grid', R_grid);
ncwrite('A408_ps_ts1_state.cdf','Z_grid', Z_grid);
ncwrite('A408_ps_ts1_state.cdf','BRRZ', BRRZ);
ncwrite('A408_ps_ts1_state.cdf','BphiRZ', BphiRZ);
ncwrite('A408_ps_ts1_state.cdf','BZRZ', BZRZ);
%}


nccreate('A408_ps_ts1_state.cdf','R_grid', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R)});
nccreate('A408_ps_ts1_state.cdf','Z_grid', 'Format', 'classic', 'Dimensions', {'dim_nz', length(EQD.Z)});
nccreate('A408_ps_ts1_state.cdf','BRRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});
nccreate('A408_ps_ts1_state.cdf','BphiRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});
nccreate('A408_ps_ts1_state.cdf','BZRZ', 'Format', 'classic', 'Dimensions', {'dim_nr', length(EQD.R),'dim_nz', length(EQD.Z)});

ncwrite('A408_ps_ts1_state.cdf','R_grid', EQD.R');
ncwrite('A408_ps_ts1_state.cdf','Z_grid', EQD.Z');
ncwrite('A408_ps_ts1_state.cdf','BRRZ', EQD.B.POL.R);
ncwrite('A408_ps_ts1_state.cdf','BphiRZ', EQD.B.TOR);
ncwrite('A408_ps_ts1_state.cdf','BZRZ', EQD.B.POL.Z);


  

