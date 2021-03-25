filename = '/home/marco/Documents/MAST/TRANSP/RUNs/29976/29976U08.CDF';
[R, Z] = nc_fluxsurfaces(filename);
GFUN = nc_read(filename, 'GFUN');
BZXR = nc_read(filename, 'BZXR');
BPHI = BZXR.*GFUN;
for k = 1:250; BPHI(k,:) = BZXR(k)*GFUN(k,:)./R(k,2:61,1); end;
 
RMAJM = nc_read(filename, 'RMAJM');
RAXIS = nc_read(filename, 'RAXIS');
BTX = nc_read(filename, 'BTX');
TIME = nc_read(filename, 'TIME');
TIME3 = nc_read(filename, 'TIME3');
BPOL = nc_read(filename, 'BPOL');
FBTX = nc_read(filename, 'FBTX');
FBPBT = nc_read(filename, 'FBPBT');
FBX = nc_read(filename, 'FBX');
BT = BTX.*FBTX;
RR = zeros(size(RMAJM));
for k = 1:length(RAXIS)
   RR(k,:) = RAXIS(k) + RMAJM(k,:) - RMAJM(k,61);
end

plot(R(100,2:61,1),BPHI(100,:))
plot(R(100,2:61,1),BPHI(100,:), RR(100,:), BT(100,:))
plot(R(100,2:61,1),BPHI(100,:), RR(100,:), BT(100,:),'+')
Q = nc_read(filename, 'Q');

addpath ('/home/marco/Documents/Fusion')
EB_kev; 70; EB_joule = energy_unit (EB_kev*1000, 1);
EB_kev= 70; EB_joule = energy_unit (EB_kev*1000, 1);
mass = 2* 1.6726E-27;
VEL = velocity (EB_joule, mass);
WBOUNCE = VEL./(0.01*RAXIS(100)*Q(100,:)).*sqrt((RR(100,62:end)-RAXIS(100))/2*RAXIS(100));


[RR, TIME, BP, BT, BTOT] = nc_b_fields_equatorial(filename);
q = 1.6022E-19;
WPREC = EB_joule./(q*0.01*RAXIS(100)^2*BP(100,:));
WCIRC = VEL/(0.01*RAXIS(100));


ND = nc_read(filename, 'ND');
mo = 4*pi*1E-7;
VALF = BT(100,62:end)./sqrt(mo*1E6*ND(100,:)*mass);

WALF = VALF./(2*0.01*RAXIS(100)*Q(100,:))


PTOWB = nc_read(filename, 'PTOWB');
PMHD_IN = nc_read(filename, 'PMHD_IN');
PPLAS = nc_read(filename, 'PPLAS');
addpath ('/home/marco/Documents/Octave')
[DR, DP] = derivative(PPLAS(100,:),0.01*(RR(100,62:end)-RAXIS(100)));
WDIAMG = mean(DP)./(q*1E6*ND(100,:).*BT(100,62:end).*(RR(100,62:end)-RAXIS(100))*0.01)

