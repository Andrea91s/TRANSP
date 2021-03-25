ascot4GC = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A407.h5');
ascot4FO = load('/home/andrea/ascot5/TRAINING_SESSION/Ascot4_fortran/29909/RUNS/A408.h5');


rho = ascot4GC.distributions.rhoDist.abscissae.dim1;

for i=1: length(rho)-1
  RHO(i) = (rho(i) + rho(i+1))/2;
end
RHODGC = ascot4GC.distributions.rhoDist.ordinate(1,:).*104;
RHODFO = ascot4FO.distributions.rhoDist.ordinate(1,:).*102;
sum(RHODGC.*ascot4GC.distributions.rhoDist.shellVolume)
sum(RHODFO.*ascot4FO.distributions.rhoDist.shellVolume)



fast = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P17/29909P17_fi_1.cdf';
big = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/P17/29909P17.CDF';

Xrho = ncread(big, 'X');
BDENS = ncread(big, 'BDENS');
DVOL = nc_read(big, 'DVOL');   
Xrho=Xrho(:,56);
BDENS = BDENS(:,56)*1e6;
DVOL = DVOL(56,:)/1e6;
NTOT_D_NBI = ncread(fast, 'NTOT_D_NBI');

RHO=RHO';
RHODGC =RHODGC';
RHODFO = RHODFO';
figure(2000)
plot(RHO, RHODGC,'b', RHO, RHODFO,'r')
legend('GC','FO')
xlim([0,1])

