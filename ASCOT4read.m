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
addpath('/home/andrea/Documents/MASTOrbit');
[EQD] = read_eqdsk('/home/andrea/ascot5/python/a5py/a5py/preprocessing/g029909.00216',0);
%[EQD] = read_eqdsk('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_matlab/g29210.00255',0);
load('RBW.mat')

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


[Rascot, Zascot] = meshgrid(Rhisto, Zhisto);
R2D_Ao = reshape(Rascot,[],1);
Z2D_Ao = reshape(Zascot,[],1);
[inside,outside] = inpolygon(R2D_Ao,Z2D_Ao,RSURF(:,end), ZSURF(:,end));
in = inside & ! outside;
R2D_A = R2D_Ao(in);
Z2D_A = Z2D_Ao(in);
BMVOL_Ao = (2*pi*( max(diff(Z2D_Ao))* max(diff(R2D_Ao)))*R2D_Ao);
BMVOL_A = BMVOL_Ao(in);

ascot4GC = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A413.h5');
#ascot4GC = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/GC/ascot.h5');
#ascot4FO = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/FO/ascot.h5');
ascot4FO = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A404.h5');

dist4agc = squeeze(ascot4GC.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adressogc = reshape(permute(dist4agc, [4,3,2,1]), [75, 50, 6771]);
dist4adressgc = dist4adressogc(:,:,in)/4744;


dist4afo = squeeze(ascot4FO.distributions.rzPitchEdist.ordinate)./6.242e18;
dist4adressofo = reshape(permute(dist4afo, [4,3,2,1]), [75, 50, 6771]);
dist4adressfo = dist4adressofo(:,:,in)/5532;




pitch4 = ascot4FO.distributions.rzPitchEdist.abscissae.dim3;
energy4 = ascot4FO.distributions.rzPitchEdist.abscissae.dim4.*6.242e18;
for i=1: length(pitch4)-1
  p4(i) = (pitch4(i) + pitch4(i+1))/2;
endfor

for i=1: length(energy4)-1
  e4(i) = (energy4(i) + energy4(i+1))/2;
endfor





[X Y] = meshgrid (linspace(min(Rhisto), max(Rhisto), length(Rhisto)),...
 linspace(min(Zhisto), max(Zhisto), length(Zhisto)));


k = sum(squeeze(sum(dist4adressgc,1)),1)*(p4(2)-p4(1)).*(e4(2)-e4(1))/2;
h = sum(squeeze(sum(dist4adressfo,1)),1)*(p4(2)-p4(1)).*(e4(2)-e4(1))/2;
FIDDAGC = griddata(R2D_A, Z2D_A, k', X, Y, 'linear');
FIDDAFO = griddata(R2D_A, Z2D_A, h', X, Y, 'linear');

u = find(isnan(FIDDAGC) == 1);
FIDDAGC(u) = 0;
clear u
u = find(isnan(FIDDAFO) == 1);
FIDDAFO(u) = 0;
clear u



figure(100)
subplot(1,3,1)
pcolor(X,Y, FIDDAGC)
hold on
plot(EQD.LCFS_R.*100, EQD.LCFS_Z.*100, 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('ASCOT4 GC, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading flat
axis equal
hold on;
plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0 4e12]);
set (gca, 'fontsize', 14)

subplot(1,3,2)
pcolor(X,Y, FIDDAFO)
hold on
plot(EQD.LCFS_R.*100, EQD.LCFS_Z.*100, 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('ASCOT4 FO, FI density cm^{-3}', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading flat
axis equal
hold on;
plot(EQD.R_axis.*100, EQD.Z_axis.*100, '+k', 'markersize', 20, 'linewidth', 3)
colormap(RBW);
colorbar('fontsize', 14)
caxis([0 4e12]);
set (gca, 'fontsize', 14)


subplot(1,3,3)
pcolor(X,Y, FIDDAFO./FIDDAGC)
hold on
plot(EQD.LCFS_R.*100, EQD.LCFS_Z.*100, 'linewidth',2,'k')
hold on
plot(RSURF(:,69), ZSURF(:,69), 'linewidth',2,'r')
title('Relative difference', 'fontsize',14)
xlabel('R (cm)','fontsize', 14)
ylabel('Z (cm)','fontsize', 14)
shading flat
axis equal
hold on;
colormap(RBW);
colorbar('fontsize', 14)
caxis([0 1.5]);
set (gca, 'fontsize', 14)

aaa = 100.*(1-(FIDDAFO./sum(FIDDAFO(:)))./(FIDDAGC./sum(FIDDAGC(:))));

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
caxis([0 2e14]);
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
caxis([0 2e14]);
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


TOTEP_ASCOT4GC = sum(FIDA4GLOBALGC(:)).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2
TOTEP_ASCOT4FO = sum(FIDA4GLOBALFO(:)).*(p4(2)-p4(1)).*(e4(2)-e4(1))/2
TOTRZ_ASCOT4GC = 2*pi*sum(sum(FIDDAGC.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))
TOTRZ_ASCOT4FO = 2*pi*sum(sum(FIDDAFO.*X)).*(X(1,2)-X(1,1))*(Y(2,1)-Y(1,1))
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

big = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/C18/29909C18.CDF';
BDENS = ncread(big, 'BDENS');
DVOL = nc_read(big, 'DVOL');
X = ncread(big, 'X');
x = X(:,56);
fid = BDENS(:,56)
vol = DVOL(56,:);
sum(fid.*vol')

FIDTRANSP = interp1(x,fid, RHO);

figure(2)
plot(RHO, FIDRHO, 'linewidth',2, 'r', RHO, FIDRHOFO, 'linewidth',2, 'b',...
RHO, FIDTRANSP*1e4, 'linewidth',2, 'g')
xlabel('\rho poloidal','fontsize', 14)
xlim([0.01,1])
ylabel('density (m^{-3})','fontsize', 14)
h = legend('GC','FO');
set (h, "fontsize", 14);
set (gca, "fontsize", 14)


fid
