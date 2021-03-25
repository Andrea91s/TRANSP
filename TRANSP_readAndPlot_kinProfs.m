% Script for reading data from a TRANSP cdf data file, plotting it, and
% writing to ASCOT inputs

% Preliminaries
%close all;


%% Define the filename of the cdf file
fn = '/home/andrea/TRANSP/TRANSP_analysis_script/RUNS/29909/29909C18.CDF';
% Alternatively use Keeling's rerun of the same case.
% NOTE: My partially careful and partially quick check showed
% no difference, i.e. that TRANSP_99999K25_2.CDF = TRANSP_99999K25.CDF
%fn = 'TRANSP_99999K25_2.CDF';


%% Read fluxes from TRANSP cdf
% Poloidal fluxes on axis
%PSI0_TR = ncread(filename, 'PSI0_TR'); % = zeros
%PSI0_DATA = ncread(filename, 'PSI0_DATA'); % = zeros
%PSI0 = ncread(filename, 'PSI0'); = 1
% Toroidal fluxes
TRFLX = ncread(fn, 'TRFLX'); % TOROIDAL FLUX [WEBERS]
%TRFMP = ncread(filename, 'TRFMP');
%TFLUX = ncread(filename, 'TFLUX');
% Poloidal fluxes
PLFLX = ncread(fn, 'PLFLX'); % POLOIDAL FLUX [Wb/rad]
%PLFMP = ncread(filename, 'PLFMP');
%PLFLXA = ncread(filename, 'PLFLXA');


%% Read time and place abscissae from TRANSP cdf
% TIME
TIME = ncread(fn,'TIME'); % TIME [SECONDS]
%TIME3 = ncread(filename,'TIME'); % Typically = TIME
% Choose time point. According to David Keeling the data reaches a nice
% equilibrium only towards the end of a run. He says 5.7s is a good
% reference time point, for scenarios A1 (K25) and A2 (K26) at least.
i_t = 56;
% PLACE
% X-coordinate should be rho_tor_norm, with sqrt taken before
% normalization. You might want to comment out the manual calculation below
% and check. Specifically, X are bin midpoints, and XB are bin boundaries,
% omitting the first, because that is always 0. Just add it manually as the
% head, also for the data signals you want to look at, and hopefully know
% the value for at the mangetic axis.
X = ncread(fn, 'X'); % x"r/a" ctr [1]
XB = ncread(fn, 'XB'); % x"r/a" bdy [1]
% Produce rho_tor, which should be the XB variable in CDF. Plot to check.
rho_tor = [0;TRFLX(:,i_t)];
%rho_tor_norm = sqrt((rho_tor - rho_tor(1))/(rho_tor(end) - rho_tor(1))); % Unnecessary because rho_tor(1) = 0
rho_tor = sqrt(rho_tor/rho_tor(end));
figure; plot([0;XB(:,i_t)], rho_tor); hold on;
xlabel('XB0'); ylabel('\rho_{tor, TRANSP}');
% Find the mid points, which should correspond to X. Plot to check.
rho_tor_mid = interp1([0;XB(:,i_t)], rho_tor, X(:,i_t));

plot(X(:,i_t), rho_tor_mid, '--');
xlabel('X'); ylabel('\rho_{tor, TRANSP}');


%% Define poloidal flux coordinate, because ASCOT uses it
rho_pol = [0;PLFLX(:,i_t)];
rho_pol = sqrt(rho_pol/rho_pol(end));
rho_pol_mid = interp1([0;XB(:,i_t)], rho_pol, X(:,i_t));
figure; plot([0;XB(:,i_t)],rho_pol,'b','linewidth',2); %hold on; plot(X(:,i_t), rho_pol_mid, '--');
xlabel('\rho_{tor, TRANSP}'); ylabel('\rho_{pol, ASCOT}');
%legend('bin boundaries', 'bin midpoints', 'Location', 'SouthEast');

%% Read shell volumes corresponding to rho_pol
% TODO: Are we sure it's rho_pol and not rho_tor?
DVOL = ncread(fn, 'DVOL'); % ZONE VOLUME [CM**3]

%% Plot normalized shell volumes
figure;
plot(rho_pol_mid, DVOL(:,i_t)./diff(rho_pol));
xlabel('\rho_{pol, TRANSP}'); ylabel('ZONE VOLUME [CM**3]');


%% Read electrons and ions from TRANSP cdf
% Read elec and D and H temperatures and densities
TE = ncread(fn, 'TE'); % ELECTRON TEMPERATURE [EV]
NE = ncread(fn, 'NE'); % ELECTRON DENSITY [N/CM**3]
TI = ncread(fn, 'TI'); % ION TEMPERATURE [EV]
ND = ncread(fn, 'ND'); % DEUTERIUM ION DENSITY [N/CM**3]
NH = ncread(fn, 'NH'); % HYDROGEN ION DENSITY [N/CM**3]
NI = ncread(fn, 'NI'); % TOTAL ION DENSITY [N/CM**3]
% Some fusion (?) and impurity densities that might be of interest:
%DNIDT = ncread(fn, 'DNIDT'); % D/DT(TOTAL ION DENSITY) [N/CM3/SEC]
%DNEDT = ncread(fn, 'DNEDT'); % D/DT(ELECTRON DENSITY) [N/CM3/SEC]
%DNIMP = ncread(fn, 'DNIMP'); % D/DT(IMPURITY DENSITY) [N/CM3/SEC]
NIMP = ncread(fn, 'NIMP'); % TOTAL IMPURITY DENSITY [N/CM**3]


%% Plot kinetic profiles
figure(1);
plot(rho_tor_mid, NE(:,i_t),'r', 'linewidth',2); hold on;
plot(rho_pol_mid, NE(:,i_t),'b', 'linewidth',2); hold on;
legend('TRANSP TOROIDAL','ASCOT POLOIDAL')
xlabel('\rho'); ylabel('electron temperature (eV)');

figure(2);
plot(rho_tor_mid, ND(:,i_t),'r', 'linewidth',2); hold on;
plot(rho_pol_mid, ND(:,i_t),'b', 'linewidth',2); hold on;
legend('TRANSP TOROIDAL','ASCOT POLOIDAL')
xlabel('\rho'); ylabel('deuterium temperature (eV)');



return


plot(rho_pol_mid, ND(:,i_t));
plot(rho_pol_mid, NH(:,i_t));
plot(rho_pol_mid, NIMP(:,i_t));
plot(rho_pol_mid, NI(:,i_t));
plot(rho_pol_mid, ND(:,i_t)+NH(:,i_t)+NIMP(:,i_t), '--');
legend('NE', 'ND', 'NH', 'NIMP', 'NI', 'ND+NH+NIMP', 'Location', 'West');
xlabel('\rho_{pol}'); ylabel('n (cm^{-3})');
box on; grid on;

figure;
plot(rho_pol_mid, TE(:,i_t)); hold on;
plot(rho_pol_mid, TI(:,i_t), '--');
legend('TE', 'TI', 'Location', 'West');
xlabel('\rho_{pol}'); ylabel('T (eV)');
box on; grid on;


%% Calculate fractions and relative deviation from flat top
numD = sum(ND(:,i_t).*DVOL(:,i_t));
numH = sum(NH(:,i_t).*DVOL(:,i_t));
H_frac = numH/(numD+numH);
disp(['Fraction of H in D+H: ', num2str(H_frac)]);
numE = sum(NE(:,i_t).*DVOL(:,i_t));
DH_frac = (numH+numD)/numE;
disp(['Fraction of D+H of E: ', num2str(DH_frac)]);
high = NI(rho_pol_mid> 0.78 & rho_pol_mid<0.8,i_t);
low = NI(1,i_t);
dev = (high-low)/high;
disp(['Deviation: ', num2str(dev)]);

%{

%% Write ASCOT plasma
% Set number of grid points and generate abscissa
Nrad = 50;
rho = linspace(rho_pol(1),rho_pol(end),Nrad);
% Load template and fill ASCOT plasma input file
pls = read_ascot('~/Documents/data/mastu/MASTU_Juan/MASTU_Juan.plasma_1d');
pls.com1 = ['# Predictive MAST-U TRANSP data for scenario K25 (A1) at ',...
            'time = 5.7 s.'];
pls.com2 = '# Read from TRANSP CDF data file 99999K25.CDF.';
pls.shot = 99999;
pls.time = 5.7;
pls.Nrad = Nrad;
pls.Nion = 2;
pls.znum = [1 1];
pls.amass = [2 1];
pls.te = interp1(rho_pol_mid,TE(:,i_t),rho,'cubic','extrap')';
% Also unit change cm^-3 --> m^-3
pls.ne = interp1(rho_pol_mid,NE(:,i_t),rho,'cubic','extrap')'*10^6;
pls.vtor = interp1(pls.rho, pls.vtor, rho)';
pls.ti = interp1(rho_pol_mid, TI(:,i_t), rho, 'cubic', 'extrap')';
% The sum of n_D, n_H and n_imp_tot equals the TRANSP n_i_tot. Since I do
% not want impurities for now, I choose to set ion densities based on
% electron density, with a D/H ratio of 9:1:
%pls.ni = [interp1(rho_pol_norm_mid,ND(:,idx_t),rho,'cubic','extrap')', ...
%          interp1(rho_pol_norm_mid,NH(:,idx_t),rho,'cubic','extrap')'] ...
%         *10^6; % Unit change cm^-3 --> m^-3
pls.ni = [pls.ne*0.9, pls.ne*0.1];
pls.zeff = interp1(pls.rho, pls.zeff, rho)';
pls.rho = rho'; % This last, because old needed above
% Write to input.plasma_1d file. NOTE: Write is now automatically switching
% to double precision. Should I perhaps do that already manually in this
% script?
write_1d_plasma_data(pls,'MASTU_K25_5700ms.plasma_1d');


%% Plot plasma for sanity check and comparison
% TRANSP vs plasma_1d
figure;
plot(rho_pol_mid, NE(:,i_t)*10^6); hold on;
plot(rho, pls.ne, '--');
xlabel('\rho_{pol}'); ylabel('n (m^{-3})');
legend('TRANSP elec', 'elec', 'Location', 'SouthWest');

% This plasma_1d vs Juan's plasma_1d
juan_pls = read_ascot('~/Documents/data/mastu/MASTU_Juan/MASTU_Juan.plasma_1d');
figure;
plot(pls.rho, pls.ne); hold on;
plot(juan_pls.rho, juan_pls.ne);
plot(pls.rho, pls.ni(:,1));
plot(juan_pls.rho, juan_pls.ni(:,1));
plot(pls.rho, pls.ni(:,2));
plot(juan_pls.rho, juan_pls.ni(:,2));
plot(pls.rho, pls.ni(:,1)+pls.ni(:,2),'--');
plot(juan_pls.rho, juan_pls.ni(:,1)+6*juan_pls.ni(:,2),'--');
xlabel('\rho_{pol}'); ylabel('n (m^{-3})');
legend('elec', 'Juan elec', 'D', 'Juan D', 'H', 'Juan H', ...
    'total ion', 'total Juan ion', 'Location', 'West');
figure;
plot(pls.rho, pls.te); hold on;
plot(juan_pls.rho, juan_pls.te);
plot(pls.rho, pls.ti,'--');
plot(juan_pls.rho, juan_pls.ti,'--');
xlabel('\rho_{pol}'); ylabel('T (eV)');
legend('elec', 'Juan elec', 'ion', 'Juan ion');


%% Read neutrals from TRANSP cdf
% Hydrogen
DN0VH = ncread(fn, 'DN0VH'); % VOL NEUTRAL DENSITY G=H [N/CM**3]
T0VH = ncread(fn, 'T0VH'); % VOL NEUTRAL TEMP G=H [EV]
DN0WH = ncread(fn, 'DN0WH'); % WALL NEUTRAL DENS G=H [N/CM**3]
T0WH = ncread(fn, 'T0WH'); % WALL NEUTRAL TEMP G=H [EV]
%N0SRC_H = ncread(filename, 'N0SRC_H'); % (= DN0WH)
% Gas flow n_0 dens is ~20 orders smaller than wall/recycling!
%N0SGF_H = ncread(filename, 'N0SGF_H'); % gas flow neutral dens G=H
% Deuterium
DN0VD = ncread(fn, 'DN0VD'); % VOL NEUTRAL DENSITY G=D [N/CM**3]
T0VD = ncread(fn, 'T0VD'); % VOL NEUTRAL TEMP G=D [EV]
DN0WD = ncread(fn, 'DN0WD'); % WALL NEUTRAL DENS G=D [N/CM**3]
T0WD = ncread(fn, 'T0WD'); % WALL NEUTRAL TEMP G=D�[EV]
%N0SRC_D = ncread(filename, 'N0SRC_D'); % recycling neutral dens G=D (= DN0WD)
% Gas flow n_0 dens is ~20 orders smaller than wall/recycling!
%N0SGF_D = ncread(filename, 'N0SGF_D'); % gas flow neutral dens G=D
%DENS0 = ncread(filename, 'DENS0'); % Empty in K25
%T0 = ncread(filename, 'T0'); % -||-


%% Plot TRANSP neutrals
figure;
semilogy(rho_pol_mid, DN0VD(:,i_t)); hold on;
semilogy(rho_pol_mid, DN0WD(:,i_t));
semilogy(rho_pol_mid, DN0VD(:,i_t)+DN0WD(:,i_t), '--');
semilogy(rho_pol_mid, DN0VH(:,i_t));
semilogy(rho_pol_mid, DN0WH(:,i_t));
semilogy(rho_pol_mid, DN0VH(:,i_t)+DN0WH(:,i_t), '--');
semilogy(rho_pol_mid, ...
    DN0VH(:,i_t)+DN0WH(:,i_t)+DN0VD(:,i_t)+DN0WD(:,i_t), ':');
xlabel('\rho_{pol}'); ylabel('n (cm^{-3})');
legend('volumetric D', 'recycled D', 'total D', ...
    'volumetric H', 'recycled H', 'total H', ...
    'total', 'Location', 'NorthWest');
figure;
plot(rho_pol_mid, T0VD(:,i_t)); hold on;
plot(rho_pol_mid, T0WD(:,i_t));
plot(rho_pol_mid, T0VH(:,i_t), '--');
plot(rho_pol_mid, T0WH(:,i_t), '--');
xlabel('\rho_{pol}'); ylabel('T (eV)');
legend('volumetric D', 'recycled D', 'volumetric H', 'recycled H', ...
    'Location', 'SouthWest');


%% Write ASCOT neutrals with extrapolation
%ntr = read_ascot('template.neutral_1d');
ntr = read_ascot('MASTU_K25_5700ms.neutral_1d');
nPoints = 50;
ntr.comment1 = ['# Predictive MAST-U TRANSP data for scenario K25 (A1) at ',...
                'time = 5.7 s.'];
ntr.comment2 = '# Read from TRANSP CDF data file 99999K25.CDF.';
rho = linspace(rho_pol(1),rho_pol(end),nPoints);
% Change to far-reaching rho that goes beyond wall, for extrapolation
rho = [rho,(rho(end)+rho((2)-rho(1))):(rho(2)-rho(1)):1.7];
ntr.nPoints = length(rho);
ntr.rho = rho';
% Density is extrapolated constantly beyond separatrix to avoid automatic
% ASCOT4 extrapolation, which would give wildly unphysical exponential
% decay. However, I do not know what kind of extrapolation would be
% physical, so this constant one is a crude approximation.
ntr.density = interp1(rho_pol_mid, ...
    DN0VH(:,i_t)+DN0WH(:,i_t)+DN0VD(:,i_t)+DN0WD(:,i_t), ...
    rho, 'linear', 'extrap')';
ntr.density(ntr.rho > 1) = ntr.density(ntr.rho == 1)*ntr.density(ntr.rho > 1)./ntr.density(ntr.rho > 1);
% Also unit change cm^-3 --> m^-3.
ntr.density = ntr.density*10^6;
% Temperature is extrapolated beyond separatrix with exponential decay,
% omitting possible negative values that the initial extrapolation by
% interp1 to the separatrix may have yielded.
ntr.temperature = interp1(rho_pol_mid, ...
    T0WD(:,i_t), ...
    rho(rho<=1), 'linear', 'extrap')';
N_temp = length(rho(rho<=1));
for i = 1:N_temp
    if(ntr.temperature(N_temp-(i-1))>0)
    f = fit([rho(N_temp-i),rho(N_temp-(i-1))]',[ntr.temperature(N_temp-i);ntr.temperature(N_temp-(i-1))],'exp1');
    break;
    end
end
ntr.temperature = [ntr.temperature(1:N_temp-(i-1)); ...
    (f.a*exp(f.b*rho(N_temp-(i-2):end)))'];
ntr.comment3 = '# Density is extrapolated linearly and temperature exponentially (omitting negative values) beyond the separatrix.';
% Write to input.plasma_1d file.
% NOTE: write_1d_neutral_data automatically switches to double precision.
write_1d_neutral_data(ntr,'MASTU_K25_5700ms.neutral_1d');


%% Plot neutral_1d
% Plot as sanity check
figure; semilogy(ntr.rho,ntr.density);
xlabel('\rho_{pol}'); ylabel('n (m^{-3})');
legend('neutral', 'Location', 'NorthWest');
figure; plot(ntr.rho,ntr.temperature);
xlabel('\rho_{pol}'); ylabel('T (eV)'); legend('neutral');


%% Read CX stuff
% POWER
PCX = ncread(fn, 'PCX'); % CHARGE EXCHANGE LOSS [WATTS/CM3] ((MAYBE))
P0NET = ncread(fn, 'P0NET'); % NET CHARGE EXCHANGE LOSS
PBCX = ncread(fn, 'PBCX'); % THERMAL ION LOSS, FAST ION CX
PBCX_D = ncread(fn, 'PBCX_D'); % THERMAL ION LOSS, CX W/ D BEAM
% FAST NEUTRALS
SBCX0 = ncread(fn, 'SBCX0'); % FAST ION CX: NEUTRALS BORN [N/CM3/SEC] ((MAYBE))
%SBCX0_D = ncread(fn, 'SBCX0_D'); % D BEAM CX: NEUTRALS BORN [N/CM3/SEC] (BEAM = beam fast particle, I think) (= SBCX0)
N0BCXD0 = ncread(fn, 'N0BCXD0'); % CX FAST NEUTRAL DENSITY (D0) [N/CM**3] ((MAYBE))
N0BD0 = ncread(fn, 'N0BD0'); % 1.GEN FAST NEUTRAL DENSITY (D0) [N/CM**3]
% (CX?) RECAPTURE
SBXR_II = ncread(fn, 'SBXR_II'); % FAST ION RECAPTURE on th.ions [N/CM3/SEC]
SBXR_IE = ncread(fn, 'SBXR_IE'); % FAST ION RECAPTURE on electrons [N/CM3/SEC]
SBXR_IZ = ncread(fn, 'SBXR_IZ'); % FAST ION RECAPTURE on impurities [N/CM3/SEC]
SBXRB = ncread(fn, 'SBXRB'); % FAST ION CX: BEAM-BEAM RECAPTURE [N/CM3/SEC]
%SBXR_I_D = ncread(fn, 'SBXR_I_D'); % D B RECAP by ioniz: th.ions (= SBXR_II)
%SBXR_E_D = ncread(fn, 'SBXR_E_D'); % D B RECAP by ioniz: electrons (= SBXR_IE)
%SBXR_Z_D = ncread(fn, 'SBXR_Z_D'); % D B RECAP by ioniz: impurities (= SBXR_IZ)
%SBXRB_D = ncread(fn, 'SBXRB_D'); % D BEAM CX: RECAPTURE: BEAM-BEAM (= SBXRB)
SBXRD = ncread(fn, 'SBXRD'); % BEAM CX: RECAPTURE BY CX W/D+ [N/CM3/SEC]
%SBXRD_D = ncread(fn, 'SBXRD_D'); % D BEAM CX: RECAPTURE BY CX W/D+ (= SBXRD)
SBXRH = ncread(fn, 'SBXRH'); % BEAM CX: RECAPTURE BY CX W/H+ [N/CM3/SEC]
%SBXRH_D = ncread(fn, 'SBXRH_D'); % D BEAM CX: RECAPTURE BY CX W/H+ (= SBXRH)
SIHALO_H = ncread(fn, 'SIHALO_H'); % BEAM HALO RECAP ION SCE G=H
SIHALO_D = ncread(fn, 'SIHALO_D'); % BEAM HALO RECAP ION SCE G=D
SEHALO = ncread(fn, 'SEHALO'); % (e-) RECAP in HALO ION SCEs
% SINKS
RSNBX_H_D = ncread(fn, 'RSNBX_H_D'); % H_0 cx sink by D beam ions [1/sec]
RSNBI_H_D = ncread(fn, 'RSNBI_H_D'); % H_0 ii sink by D beam ions [1/sec]
RSNBX_D_D = ncread(fn, 'RSNBX_D_D'); % D_0 cx sink by D beam ions [1/sec] ((MAYBE, but what is this?))
RSNBI_D_D = ncread(fn, 'RSNBI_D_D'); % D_0 ii sink by D beam ions [1/sec]
SNCX_D = ncread(fn, 'SNCX_D'); % CX sink rate beam D [N/CM3/SEC] ((MAYBE, but is this beam ionization?))
SNCXMC_D = ncread(fn, 'SNCXMC_D'); % MC CX sink rate beam D,orbit [N/CM3/SEC]
% GAS FLOW ((Recycling data better. Check corresponding neutral densities!))
% SFCXSGF = ncread(fn, 'SFCXSGF'); % gas flow (e-)=> FAST ION CX [N/CM3/SEC] (Sum H+D below)
% T0CX_GFH = ncread(fn, 'T0CX_GFH'); % CX NEUTRAL TEMP. gas flow H [EV]
% SFCXGF_H = ncread(fn, 'SFCXGF_H'); % H gas (e-)=> FAST ION CX [N/CM3/SEC]
% T0CX_GFD = ncread(fn, 'T0CX_GFD'); % CX NEUTRAL TEMP. gas flow D [EV]
% SFCXGF_D = ncread(fn, 'SFCXGF_D'); % D gas (e-)=> FAST ION CX [N/CM3/SEC]
% RECYCLING
SFCXSRC = ncread(fn, 'SFCXSRC'); % recycling (e-)=> FAST ION CX [N/CM3/SEC] (Sum H+D below) ((MAYBE))
T0CX_RCH = ncread(fn, 'T0CX_RCH'); % CX NEUTRAL TEMP. recyc. H [EV]
SFCXRC_H = ncread(fn, 'SFCXRC_H'); % H recyc (e-)=> FAST ION CX [N/CM3/SEC]
T0CX_RCD = ncread(fn, 'T0CX_RCD'); % CX NEUTRAL TEMP. recyc. D [EV]
SFCXRC_D = ncread(fn, 'SFCXRC_D'); % D recyc (e-)=> FAST ION CX [N/CM3/SEC]
% BEAMS
PCX_HALO = ncread(fn, 'PCX_HALO'); % beam halo driven cx power [WATTS/CM3]
SFCXHALO = ncread(fn, 'SFCXHALO'); % HALO NEUTRALS (e-)=> FAST ION CX [N/CM3/SEC]
% ONLY TIME DEP STUFF
%BPCXI = ncread(fn, 'BPCXI'); % FAST ION POWER TO CX (INT) (= BPCXI_D below)
%BPCXX = ncread(fn, 'BPCXX'); % FAST ION POWER TO CX (EXT) (= BPCXX_D below)
BSNXI = ncread(fn, 'BSNXI'); % FAST ION CX SINK (INT)
BSNXO = ncread(fn, 'BSNXO'); % FAST ION CX SINK (EXT)
SBCXX = ncread(fn, 'SBCXX'); % CX FAST ION LOSS [N/SEC]
%SBCXX_D = ncread(fn, 'SBCXX_D'); % CX D BEAM ION LOSS [N/SEC] (= SBCXX)
NCX0_D = ncread(fn, 'NCX0_D'); % # CX events D orbiting [N]
SFRCAP = ncread(fn, 'SFRCAP'); % FAST ION CX RECAPTURE [N/SEC] (likely complemented by EI, II and ZI ionization)
P0ESC = ncread(fn, 'P0ESC'); % NEUTRAL POWER ESCAPED [WATTS] (not fast ion -dedicated)
BPHTO = ncread(fn, 'BPHTO'); % TOTAL FAST ION HEATING [WATTS]
BPTI = ncread(fn, 'BPTI'); % BEAM POWER TO IONS [WATTS] (this...)
BPTE = ncread(fn, 'BPTE'); % BEAM POWER TO ELECTRONS [WATTS] (...and this perhaps complemented by P to impurities to give BPHTO)
%BPTI_D = ncread(fn, 'BPTI_D'); % D BEAM POWER TO IONS [WATTS] (= BPTI)
%BPTE_D = ncread(fn, 'BPTE_D'); % D BEAM POWER TO ELECTRONS [WATTS] (= BPTE)
BPSHI = ncread(fn, 'BPSHI'); % FAST ION SHINE-THRU POWER [WATTS]
BPCAP = ncread(fn, 'BPCAP'); % BEAM POWER CAPTURED [WATTS]
BPCXI_D = ncread(fn, 'BPCXI_D'); % D BEAM POWER TO CX (INT) [WATTS]
BPCXX_D = ncread(fn, 'BPCXX_D'); % D BEAM POWER TO CX (EXT) [WATTS]
BPCI0_D = ncread(fn, 'BPCI0_D'); % D BEAM CX SCE POWER (INT) [WATTS]
BPCX0_D = ncread(fn, 'BPCX0_D'); % D BEAM CX SCE POWER (EXT) [WATTS]
BPCRI_D = ncread(fn, 'BPCRI_D'); % D BEAM CX RECAPTURE (INT) [WATTS]
BPCRX_D = ncread(fn, 'BPCRX_D'); % D BEAM CX RECAPTURE (EXT) [WATTS]
%PBSHINE_D = ncread(fn, 'PBSHINE_D'); % DBEAM SHINE-THRU POWER [WATTS] (= BPSHI)


%% Plot TRANSP CX stuff for benchmark
% POWER
% figure;
% semilogy(rho_pol_norm_mid,PCX(:,i_t)); hold on;
% xlabel('\rho_{pol}'); ylabel('P (Wcm^{-3})');
% legend('PCX = CHARGE EXCHANGE LOSS');

% FAST NEUTRALS BORN
figure;
semilogy(rho_pol_mid,SBCX0(:,i_t)); hold on;
semilogy(rho_pol_mid,SFCXRC_H(:,i_t));
semilogy(rho_pol_mid,SFCXRC_D(:,i_t));
semilogy(rho_pol_mid,SFCXSRC(:,i_t), '--');
semilogy(rho_pol_mid,SBCX0(:,i_t)+SFCXSRC(:,i_t), '-.');
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('FAST ION CX: NEUTRALS BORN', ...
       'H recyc (e-)=> FAST ION CX', ...
       'D recyc (e-)=> FAST ION CX', ...
       'recycling (e-)=> FAST ION CX', ...
       'sum', 'Location', 'SouthEast');
% (CX?) RECAPTURE
figure;
plot(rho_pol_mid,SBXR_II(:,i_t)); hold on;
plot(rho_pol_mid,SBXR_IE(:,i_t));
plot(rho_pol_mid,SBXR_IZ(:,i_t));
plot(rho_pol_mid,SBXRB(:,i_t));
plot(rho_pol_mid,SBXRD(:,i_t));
plot(rho_pol_mid,SBXRH(:,i_t));
plot(rho_pol_mid,SBXR_II(:,i_t)+SBXR_IE(:,i_t)+SBXR_IZ(:,i_t)+SBXRB(:,i_t)+SBXRD(:,i_t)+SBXRH(:,i_t), '--');
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('FAST ION RECAPTURE on th.ions', ...
       'FAST ION RECAPTURE on electrons', ...
       'FAST ION RECAPTURE on impurities', ...
       'FAST ION CX: BEAM-BEAM RECAPTURE', ...
       'BEAM CX: RECAPTURE BY CX W/D+', ...
       'BEAM CX: RECAPTURE BY CX W/H+', ...
       'sum');
% FAST IONS BORN
figure;
semilogy(rho_pol_mid,SBXR_II(:,i_t)+SBXR_IE(:,i_t)+SBXR_IZ(:,i_t)+...
                          SBXRB(:,i_t)+SBXRD(:,i_t)+SBXRH(:,i_t)); hold on;
semilogy(rho_pol_mid,SNCX_D(:,i_t));
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('TRANSP: all RECAPTURE', ...
       'TRANSP: SNCX_D = CX sink rate beam D', ...
       'Location', 'NorthWest');
% 0D QUANTITIES (ONLY TIME DEP)
disp(['CX FAST ION LOSS = ', num2str(SBCXX(i_t)), ' 1/s']);
disp(['NEUTRAL POWER ESCAPED = ', num2str(P0ESC(i_t)), ' W']);
disp(['BEAM POWER TO IONS: ', num2str(BPTI(i_t)), ' WATTS']);
disp(['BEAM POWER TO ELECTRONS: ', num2str(BPTE(i_t)), ' WATTS']);
disp(['Sum of above two: ', num2str(BPTI(i_t)+BPTE(i_t)), ' WATTS']);
disp(['TOTAL FAST ION HEATING: ', num2str(BPHTO(i_t)), ' WATTS']);
disp(['FAST ION SHINE-THRU POWER = ', num2str(BPSHI(i_t)), ' W']);
disp(['BEAM POWER CAPTURED = ', num2str(BPCAP(i_t)), ' W']);
disp(['D BEAM POWER TO CX (INT + EXT) = ', num2str(BPCXI_D(i_t)+BPCXX_D(i_t)), ' W']);
disp(['D BEAM CX SCE POWER (INT + EXT) = ', num2str(BPCI0_D(i_t)+BPCX0_D(i_t)), ' W']);
disp(['D BEAM CX RECAPTURE (INT + EXT) = ', num2str(BPCRI_D(i_t)+BPCRX_D(i_t)), ' W']);
disp(['Comparison of above 3, 1st vs (2nd - 3rd) = ', ...
      num2str(BPCXI_D(i_t)+BPCXX_D(i_t)), ...
      ' vs ', ...
      num2str((BPCI0_D(i_t)+BPCX0_D(i_t))- ...
              (BPCRI_D(i_t)+BPCRX_D(i_t))), ' W']);
disp(['Time-averaged comparison of above 3, 1st vs (2nd - 3rd) = ', ...
      num2str(mean(BPCXI_D(i_t-2:i_t+2)+BPCXX_D(i_t-2:i_t+2))), ...
      ' vs ', ...
      num2str(mean((BPCI0_D(i_t-2:i_t+2)+BPCX0_D(i_t-2:i_t+2))- ...
                   (BPCRI_D(i_t-2:i_t+2)+BPCRX_D(i_t-2:i_t+2)))), ' W']);
disp('Total recapture rates (s^{-1})');
disp(['FAST ION RECAPTURE on th.ions: ', num2str(sum(SBXR_II(:,i_t).*DVOL(:,i_t)))]);
disp(['FAST ION RECAPTURE on electrons: ', num2str(sum(SBXR_IE(:,i_t).*DVOL(:,i_t)))]);
disp(['FAST ION RECAPTURE on impurities: ', num2str(sum(SBXR_IZ(:,i_t).*DVOL(:,i_t)))]);
disp(['FAST ION CX: BEAM-BEAM RECAPTURE: ', num2str(sum(SBXRB(:,i_t).*DVOL(:,i_t)))]);
disp(['BEAM CX: RECAPTURE BY CX W/D+: ', num2str(sum(SBXRD(:,i_t).*DVOL(:,i_t)))]);
disp(['BEAM CX: RECAPTURE BY CX W/H+: ', num2str(sum(SBXRH(:,i_t).*DVOL(:,i_t)))]);
disp(['Fraction of recaptures on impurities: ', num2str((sum(SBXR_IZ(:,i_t).*DVOL(:,i_t)))/ ...
    (sum(SBXR_II(:,i_t).*DVOL(:,i_t)) + sum(SBXR_IE(:,i_t).*DVOL(:,i_t)) + ...
     sum(SBXR_IZ(:,i_t).*DVOL(:,i_t)) + sum(SBXRB(:,i_t).*DVOL(:,i_t)) + ...
     sum(SBXRD(:,i_t).*DVOL(:,i_t)) + sum(SBXRH(:,i_t).*DVOL(:,i_t))))]);
disp(['Fraction of recaptures on electrons: ', num2str((sum(SBXR_IE(:,i_t).*DVOL(:,i_t)))/ ...
    (sum(SBXR_II(:,i_t).*DVOL(:,i_t)) + sum(SBXR_IE(:,i_t).*DVOL(:,i_t)) + ...
     sum(SBXR_IZ(:,i_t).*DVOL(:,i_t)) + sum(SBXRB(:,i_t).*DVOL(:,i_t)) + ...
     sum(SBXRD(:,i_t).*DVOL(:,i_t)) + sum(SBXRH(:,i_t).*DVOL(:,i_t))))]);



%% Plot more TRANSP CX stuff
% POWER
figure;
semilogy(rho_pol_mid,PCX(:,i_t)); hold on;
semilogy(rho_pol_mid,P0NET(:,i_t), '--');
semilogy(rho_pol_mid,PBCX(:,i_t));
semilogy(rho_pol_mid,PBCX_D(:,i_t), '--');
xlabel('\rho_{pol}'); ylabel('P (Wcm^{-3})');
legend('PCX = CHARGE EXCHANGE LOSS', 'net CX loss', ...
    'th ion loss by fast ion CX', 'th ion loss by beam D CX');
% FAST NEUTRALS
figure;
semilogy(rho_pol_mid,SBCX0(:,i_t)); hold on;
xlabel('\rho_{pol}'); ylabel('rate (cm^{-3}s^{-1})');
legend('SBCX0 = FAST ION CX: NEUTRALS BORN');
figure;
semilogy(rho_pol_mid,N0BCXD0(:,i_t)); hold on;
semilogy(rho_pol_mid,N0BD0(:,i_t), '--');
xlabel('\rho_{pol}');  ylabel('density (cm^{-3})');
legend('N0BCXD0 = CX FAST NEUTRAL DENSITY', '1. gen fast ntr');
% SINKS
figure;
semilogy(rho_pol_mid,RSNBX_H_D(:,i_t)); hold on;
semilogy(rho_pol_mid,RSNBI_H_D(:,i_t));
semilogy(rho_pol_mid,RSNBX_D_D(:,i_t), '--');
semilogy(rho_pol_mid,RSNBI_D_D(:,i_t), '--');
xlabel('\rho_{pol}');  ylabel('rate (s^{-1})');
legend('H_0 cx sink by D beam ions', 'H_0 ii sink by D beam ions', ...
       'D_0 cx sink by D beam ions', 'D_0 ii sink by D beam ions');
figure;
semilogy(rho_pol_mid,SNCX_D(:,i_t)); hold on;
semilogy(rho_pol_mid,SNCXMC_D(:,i_t), '--');
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('CX sink rate beam D', 'MC CX beam D');
% GAS FLOW ((Recycling data better. Check corresponding neutral densities!))
%figure;
%semilogy(rho_pol_norm_mid,SFCXGF_D(:,i_t)); hold on;
%semilogy(rho_pol_norm_mid,SFCXGF_H(:,i_t));
%semilogy(rho_pol_norm_mid,SFCXSGF(:,i_t), '--');
% Check that total is correct
%semilogy(rho_pol_norm_mid,SFCXGF_D(:,i_t)+SFCXGF_H(:,i_t));
%xlabel('\rho_{pol}');  ylabel('gas flow - fast ion CX rate (cm^{-3}s^{-1})');
%legend('D gas', 'H gas', 'total');
%figure;
%plot(rho_pol_norm_mid,T0CX_GFD(:,i_t)); hold on;
%plot(rho_pol_norm_mid,T0CX_GFH(:,i_t), '--');
%xlabel('\rho_{pol}');  ylabel('ntr gas temp (eV)');
%legend('D', 'H');
% RECYCLING
figure;
semilogy(rho_pol_mid,SFCXRC_D(:,i_t)); hold on;
semilogy(rho_pol_mid,SFCXRC_H(:,i_t));
semilogy(rho_pol_mid,SFCXSRC(:,i_t), '--');
% Check that total is correct
%semilogy(rho_pol_norm_mid,SFCXRC_D(:,i_t)+SFCXRC_H(:,i_t));
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('SFCXRC_D = H recyc (e-)=> FAST ION CX', ...
       'SFCXRC_H = D recyc (e-)=> FAST ION CX', ...
       'SFCXSRC = recycling (e-)=> FAST ION CX');
figure;
plot(rho_pol_mid,T0CX_RCD(:,i_t)); hold on;
plot(rho_pol_mid,T0CX_RCH(:,i_t), '--');
xlabel('\rho_{pol}'); ylabel('ntr recycling temp (eV)');
legend('D', 'H');
% BEAMS
% figure;
% plot(rho_pol_norm_mid,SFCXHALO(:,i_t)); hold on;
% xlabel('\rho_{pol}'); ylabel('rate (cm^{-3}s^{-1})');
% legend('halo ntr - fast ion CX');
% figure;
% plot(rho_pol_norm_mid,PCX_HALO(:,i_t)); hold on;
% xlabel('\rho_{pol}'); ylabel('P (Wcm^{-3})');
% legend('Beam halo CX');
% ONLY TIME DEP STUFF
% Nothing yet


%% Read other relevant quantities from TRANSP CDF
% BEAM DEPOSITION
SDBBI = ncread(fn, 'SDBBI'); % BEAM DEPOSITION: BEAM-BEAM II [N/CM3/SEC]
SDBBX = ncread(fn, 'SDBBX'); % BEAM DEPOSITION: BEAM-BEAM CX [N/CM3/SEC]
SDB_II = ncread(fn, 'SDB_II'); % BEAM DEP: ioniz. on therm. ions [N/CM3/SEC]
SDB_IE = ncread(fn, 'SDB_IE'); % BEAM DEP: ioniz. on electrons [N/CM3/SEC]
SDB_IZ = ncread(fn, 'SDB_IZ'); % BEAM DEP: ioniz. on impurities [N/CM3/SEC]
BDEP_D = ncread(fn, 'BDEP_D'); % D BEAM DEPOSITION (TOTAL) [N/CM3/SEC]
SDEP_D = ncread(fn, 'SDEP_D'); % D BEAM ORBIT AV DEP (TOTAL) [N/CM3/SEC] (orbit averaged?)
%SDBBI_D = ncread(fn, 'SDBBI_D'); % D BEAM DEPOSITION: BEAM-BEAM II [N/CM3/SEC] (= SDBBI)
%SDBBX_D = ncread(fn, 'SDBBX_D'); % D BEAM DEPOSITION: BEAM-BEAM CX [N/CM3/SEC] (= SDBBX)
%SDBII_D = ncread(fn, 'SDBII_D'); % D BEAM DEP: IONIZ. on therm.ions [N/CM3/SEC] (= SDB_II)
%SDBIE_D = ncread(fn, 'SDBIE_D'); % D BEAM DEP: IONIZ. on electrons [N/CM3/SEC] (= SDB_IE)
%SDBIZ_D = ncread(fn, 'SDBIZ_D'); % D BEAM DEP: IONIZ. on impurities [N/CM3/SEC] (= BDEP_D)
SDCXD = ncread(fn, 'SDCXD'); % BEAM DEPOSITION: CX W/D+ IONS [N/CM3/SEC]
%SDCXD_D = ncread(fn, 'SDCXD_D'); % D BEAM DEPOSITION: CX W/D+ IONS [N/CM3/SEC] (= SDCXD)
SDCXH = ncread(fn, 'SDCXH'); % BEAM DEPOSITION: CX W/H+ IONS [N/CM3/SEC]
%SDCXH_D = ncread(fn, 'SDCXH_D'); % D BEAM DEPOSITION: CX W/H+ IONS [N/CM3/SEC] (= SDCXH)
% BEAM DENSITY
%BDENS2_D = ncread(fn, 'BDENS2_D'); % D Beam ion density, GC
NB_F1_D = ncread(fn, 'NB_F1_D'); % density: full energy D beam
NB_F2_D = ncread(fn, 'NB_F2_D'); % density: half energy D beam
NB_F3_D = ncread(fn, 'NB_F3_D'); % density: 1/3 energy D beam
BDENS = ncread(fn, 'BDENS'); % BEAM ION DENSITY [N/CM**3] (NOT exactly sum of above! What could be missing?)
%BDENS_D = ncread(fn, 'BDENS_D'); % D BEAM ION DENSITY (= BDENS)
%MCDENS_D = ncread(fn, 'MCDENS_D'); % D BEAM ION DENSITY (MC LIST)
% OTHER
Q = ncread(fn, 'Q'); % Q PROFILE
ZEFFP = ncread(fn, 'ZEFFP'); % PLASMA COMPOSITION ZEFF PROFILE
BTRAP_D = ncread(fn, 'BTRAP_D'); % D beam ions banana fraction
SBTH = ncread(fn, 'SBTH'); % FAST ION THERMALIZATION SOURCE
SFETO = ncread(fn, 'SFETO'); % ELECTRONS -> FAST NEUTRALS
S0VLE = ncread(fn, 'S0VLE'); % TOTAL NEUTRAL VOL SCE
Q0 = ncread(fn, 'Q0'); % Q ON AXIS
PINJ = ncread(fn, 'PINJ'); % BEAM POWER INJECTED [WATTS]


%% Plot other relevant quantities from TRANSP cdf
% BEAM DEPOSITION
figure;
plot(rho_pol_mid,SDBBI(:,i_t)); hold on;
plot(rho_pol_mid,SDBBX(:,i_t));
plot(rho_pol_mid,SDB_II(:,i_t));
plot(rho_pol_mid,SDB_IE(:,i_t));
plot(rho_pol_mid,SDB_IZ(:,i_t));
plot(rho_pol_mid,SDCXD(:,i_t));
plot(rho_pol_mid,SDCXH(:,i_t));
plot(rho_pol_mid,SDBBI(:,i_t)+SDBBX(:,i_t)+SDB_II(:,i_t)+SDB_IE(:,i_t)+SDB_IZ(:,i_t)+SDCXD(:,i_t)+SDCXH(:,i_t));
plot(rho_pol_mid,BDEP_D(:,i_t), '--');
plot(rho_pol_mid,SDEP_D(:,i_t), ':');
xlabel('\rho_{pol}');  ylabel('rate (cm^{-3}s^{-1})');
legend('BEAM DEPOSITION: BEAM-BEAM II', ...
       'BEAM DEPOSITION: BEAM-BEAM CX', ...
       'BEAM DEP: ioniz. on therm. ions', ...
       'BEAM DEP: ioniz. on electrons', ...
       'BEAM DEP: ioniz. on impurities', ...
       'BEAM DEPOSITION: CX W/D+ IONS', ...
       'BEAM DEPOSITION: CX W/H+ IONS', ...
       'sum of above', ...
       'D BEAM DEPOSITION (TOTAL)', ...
       'D BEAM ORBIT AV DEP (TOTAL)');
% FRACTION DEPOSITED ON ELECTRONS
disp(['Fraction of deposition on electrons: ', num2str((sum(SDB_IE(:,i_t).*DVOL(:,i_t)))/ ...
    (sum(SDBBI(:,i_t).*DVOL(:,i_t)) + sum(SDBBX(:,i_t).*DVOL(:,i_t)) + ...
     sum(SDB_II(:,i_t).*DVOL(:,i_t)) + sum(SDB_IE(:,i_t).*DVOL(:,i_t)) + sum(SDB_IZ(:,i_t).*DVOL(:,i_t)) + ...
     sum(SDCXD(:,i_t).*DVOL(:,i_t)) + sum(SDCXH(:,i_t).*DVOL(:,i_t))))]);
% BEAM DENSITY
figure;
plot(rho_pol_mid, NB_F1_D(:,i_t), ':'); hold on;
plot(rho_pol_mid, NB_F2_D(:,i_t), ':');
plot(rho_pol_mid, NB_F3_D(:,i_t), ':');
plot(rho_pol_mid, NB_F1_D(:,i_t)+NB_F2_D(:,i_t)+NB_F3_D(:,i_t));
plot(rho_pol_mid, BDENS(:,i_t));
xlabel('\rho_{pol}');  ylabel('density (cm^{-3})');
legend('density: full energy D beam', ...
       'density: half energy D beam', ...
       'density: 1/3 energy D beam', ...
       'sum of above', ...
       'BEAM ION DENSITY');


%% Read beam deposition, from dedicated cdf variables
% Read beam deposition for both beams
BDEP01_TOT = ncread('TRANSP_99999K25.CDF', 'BDEP01_TOT'); % bdep: Beam#01(D),total depositio [N/CM3/SEC]
BDEP01_E1 = ncread('TRANSP_99999K25.CDF', 'BDEP01_E1'); % bdep: Beam#01(D), E-frac no.1 [N/CM3/SEC]
BDEP01_E2 = ncread('TRANSP_99999K25.CDF', 'BDEP01_E2'); % bdep: Beam#01(D), E-frac no.2 [N/CM3/SEC]
BDEP01_E3 = ncread('TRANSP_99999K25.CDF', 'BDEP01_E3');% bdep: Beam#01(D), E-frac no.3 [N/CM3/SEC]
PBTOT01 = ncread('TRANSP_99999K25.CDF', 'PBTOT01');% Beam#01(D), total power [WATTS/CM3]
BDEP02_TOT = ncread('TRANSP_99999K25.CDF', 'BDEP02_TOT'); % nb: Beam#02(D), total density [N/CM3/SEC]
BDEP02_E1 = ncread('TRANSP_99999K25.CDF', 'BDEP02_E1'); % nb: Beam#02(D), E-frac no.1 [N/CM3/SEC]
BDEP02_E2 = ncread('TRANSP_99999K25.CDF', 'BDEP02_E2'); % nb: Beam#02(D), E-frac no.2 [N/CM3/SEC]
BDEP02_E3 = ncread('TRANSP_99999K25.CDF', 'BDEP02_E3'); % nb: Beam#02(D), E-frac no.3 [N/CM3/SEC]
PBTOT02 = ncread('TRANSP_99999K25.CDF', 'PBTOT02');% Beam#02(D), total power [WATTS/CM3]


%% Calc E fracs from beam deposition by integration over shell volumes
% Q: Why is BDEP01_TOT + BDEP02_TOT != BDEP_D?
q_e = 1.602176565e-19;
prt_fracs_B1 = [sum(BDEP01_E1(:,i_t).*DVOL(:,i_t)), ...
                sum(BDEP01_E2(:,i_t).*DVOL(:,i_t)), ...
                sum(BDEP01_E3(:,i_t).*DVOL(:,i_t))];
prt_sum_B1 = sum(prt_fracs_B1);
prt_tot_B1 = sum(BDEP01_TOT(:,i_t).*DVOL(:,i_t));
P_fracs_B1 = prt_fracs_B1.*[75*q_e,(75/2)*q_e,(75/3)*q_e];
P_sum_B1 = sum(P_fracs_B1);
prt_fracs_B2 = [sum(BDEP02_E1(:,i_t).*DVOL(:,i_t)), ...
                sum(BDEP02_E2(:,i_t).*DVOL(:,i_t)), ...
                sum(BDEP02_E3(:,i_t).*DVOL(:,i_t))];
prt_sum_B2 = sum(prt_fracs_B2);
prt_tot_B2 = sum(BDEP02_TOT(:,i_t).*DVOL(:,i_t));
P_fracs_B2 = prt_fracs_B2.*[75*q_e,(75/2)*q_e,(75/3)*q_e];
P_sum_B2 = sum(P_fracs_B2);
disp(['Sum of parts vs total, beam 1: ', num2str(prt_sum_B1), ' vs ', num2str(prt_tot_B1)]);
disp(['Sum of parts vs total, beam 2: ', num2str(prt_sum_B2), ' vs ', num2str(prt_tot_B2)]);
disp(['P fractions of beam 1: ', num2str(P_fracs_B1/P_sum_B1)]);
disp(['P fractions of beam 2: ', num2str(P_fracs_B2/P_sum_B2)]);


%% Plot beam deposition, from dedicated cdf variables
% Check that sum and total agree and compare beam 1 and 2
figure;
plot(rho_pol_mid, BDEP01_TOT(:,i_t)); hold on;
%plot(rho_pol_norm_mid, BDEP01_E1(:,i_t)+BDEP01_E2(:,i_t)+BDEP01_E3(:,i_t), '--'); % (~= BDEP01_TOT)
plot(rho_pol_mid, BDEP02_TOT(:,i_t));
%plot(rho_pol_norm_mid, BDEP02_E1(:,i_t)+BDEP02_E2(:,i_t)+BDEP02_E3(:,i_t), '--'); % (~= BDEP02_TOT)
plot(rho_pol_mid, BDEP01_TOT(:,i_t)+BDEP02_TOT(:,i_t), '--');
legend('BDEP01_TOT = bdep: Beam#01(D),total depositio', ...
       'BDEP02_TOT = nb: Beam#02(D), total density', ...
       'sum');
xlabel('rho'); ylabel('rate (cm^{-3}s^{-1})');
% Plot powers
figure;
plot(rho_pol_mid, PBTOT01(:,i_t)); hold on;
plot(rho_pol_mid, PBTOT02(:,i_t));
plot(rho_pol_mid, PBTOT01(:,i_t)+PBTOT02(:,i_t), '--');
legend('PBTOT01 = Beam#01(D), total power', ...
       'PBTOT02 = Beam#02(D), total power', ...
       'sum');
xlabel('rho'); ylabel('P (Wcm^{-3})');


%% Read beam density, from dedicated cdf variables
% Read beam density for both beams
NB01_TOT = ncread('TRANSP_99999K25.CDF', 'NB01_TOT'); % nb: Beam#01(D), total density [N/CM**3]
NB01_E1 = ncread('TRANSP_99999K25.CDF', 'NB01_E1'); % nb: Beam#01(D), E-frac no.1 [N/CM**3]
NB01_E2 = ncread('TRANSP_99999K25.CDF', 'NB01_E2'); % nb: Beam#01(D), E-frac no.2 [N/CM**3]
NB01_E3 = ncread('TRANSP_99999K25.CDF', 'NB01_E3'); % nb: Beam#01(D), E-frac no.3 [N/CM**3]
NB02_TOT = ncread('TRANSP_99999K25.CDF', 'NB02_TOT'); % nb: Beam#02(D), total density [N/CM**3]
NB02_E1 = ncread('TRANSP_99999K25.CDF', 'NB02_E1'); % nb: Beam#02(D), E-frac no.1 [N/CM**3]
NB02_E2 = ncread('TRANSP_99999K25.CDF', 'NB02_E2'); % nb: Beam#02(D), E-frac no.2 [N/CM**3]
NB02_E3 = ncread('TRANSP_99999K25.CDF', 'NB02_E3'); % nb: Beam#02(D), E-frac no.3 [N/CM**3]


%% Plot beam density, from dedicated cdf variables
% Calculate energy fractions by averaging
% WRONG: I need to integrate over rho volumes and calculate the
% fractions from the total density of each fraction!
E_fracs_B1 = zeros(1,3); % WRONG
E_fracs_B1(1) = mean(NB01_E1(:,i_t)./(NB01_E1(:,i_t)+NB01_E2(:,i_t)+NB01_E3(:,i_t))); % WRONG
E_fracs_B1(2) = mean(NB01_E2(:,i_t)./(NB01_E1(:,i_t)+NB01_E2(:,i_t)+NB01_E3(:,i_t))); % WRONG
E_fracs_B1(3) = mean(NB01_E3(:,i_t)./(NB01_E1(:,i_t)+NB01_E2(:,i_t)+NB01_E3(:,i_t))); % WRONG
E_fracs_B2 = zeros(1,3); % WRONG
E_fracs_B2(1) = mean(NB02_E1(:,i_t)./(NB02_E1(:,i_t)+NB02_E2(:,i_t)+NB02_E3(:,i_t))); % WRONG
E_fracs_B2(2) = mean(NB02_E2(:,i_t)./(NB02_E1(:,i_t)+NB02_E2(:,i_t)+NB02_E3(:,i_t))); % WRONG
E_fracs_B2(3) = mean(NB02_E3(:,i_t)./(NB02_E1(:,i_t)+NB02_E2(:,i_t)+NB02_E3(:,i_t))); % WRONG
disp(['Energy beam density fractions of beam 1: ', num2str(E_fracs_B1)]); % WRONG
disp(['Energy beam density fractions of beam 2: ', num2str(E_fracs_B2)]); % WRONG
% Check that sum and total agree and compare beam 1 and 2
figure;
plot(rho_pol_mid, NB01_TOT(:,i_t)); hold on;
%plot(rho_pol_norm_mid, NB01_E1(:,i_t)+NB01_E2(:,i_t)+NB01_E3(:,i_t), '--'); % (= NB01_TOT)
plot(rho_pol_mid, NB02_TOT(:,i_t));
%plot(rho_pol_norm_mid, NB02_E1(:,i_t)+NB02_E2(:,i_t)+NB02_E3(:,i_t), '--'); % (= NB02_TOT)
plot(rho_pol_mid, NB01_TOT(:,i_t)+NB02_TOT(:,i_t), '--');
legend('NB01_TOT = nb: Beam#01(D), total density', ...
       'NB02_TOT = nb: Beam#02(D), total density', ...
       'sum');
xlabel('rho'); ylabel('density (cm^{-3})');


%% Read yet another set of energy-fraction-specific beam density variables and calc energy fractions
NB_F1_D = ncread('TRANSP_99999K25.CDF', 'NB_F1_D'); % density: full energy D beam [N/CM**3]
NB_F2_D = ncread('TRANSP_99999K25.CDF', 'NB_F2_D'); % density: half energy D beam [N/CM**3]
NB_F3_D = ncread('TRANSP_99999K25.CDF', 'NB_F3_D'); % density: 1/3 energy D beam [N/CM**3]
% WRONG: I need to integrate over rho volumes and calculate the
% fractions from the total density of each fraction!
E_fracs = zeros(1,3); % WRONG
E_fracs(1) = mean(NB_F1_D(:,i_t)./(NB_F1_D(:,i_t)+NB_F2_D(:,i_t)+NB_F3_D(:,i_t))); % WRONG
E_fracs(2) = mean(NB_F2_D(:,i_t)./(NB_F1_D(:,i_t)+NB_F2_D(:,i_t)+NB_F3_D(:,i_t))); % WRONG
E_fracs(3) = mean(NB_F3_D(:,i_t)./(NB_F1_D(:,i_t)+NB_F2_D(:,i_t)+NB_F3_D(:,i_t))); % WRONG
disp(['Energy fractions of beams: ', num2str(E_fracs)]); % WRONG


%% Read plasma current
CUR = ncread('TRANSP_99999K25.CDF', 'CUR'); % TOTAL PLASMA CURRENT [AMPS/CM2]
CUROH = ncread('TRANSP_99999K25.CDF', 'CUROH'); % OHMIC PLASMA CURRENT [AMPS/CM2]
CURB = ncread('TRANSP_99999K25.CDF', 'CURB'); % BEAM DRIVEN CURRENT [AMPS/CM2]
PCUR = ncread('TRANSP_99999K25.CDF', 'PCUR'); % MEASURED PLASMA CURRENT [AMPS]
PCUREQ = ncread('TRANSP_99999K25.CDF', 'PCUREQ'); % EQ PLASMA CURRENT [AMPS]
PCURC = ncread('TRANSP_99999K25.CDF', 'PCURC'); % CALCULATED PLASMA CURRENT [AMPS]


%% Plot and print plasma current
figure;
plot(rho_pol_mid, CUR(:,i_t)); hold on;
plot(rho_pol_mid, CUROH(:,i_t));
plot(rho_pol_mid, CURB(:,i_t));
plot(rho_pol_mid, CUROH(:,i_t)+CURB(:,i_t), '--');
legend('TOTAL PLASMA CURRENT', 'OHMIC PLASMA CURRENT', 'BEAM DRIVEN CURRENT', 'sum of two previous');
xlabel('\rho_{pol}'); ylabel('current density (Acm^{-2})');
disp(['Plasma current in A, measured, equilibrium and calculated: ', ...
    num2str(PCUR(i_t)), ', ', num2str(PCUREQ(i_t)), ', ', num2str(PCURC(i_t))]);


%}