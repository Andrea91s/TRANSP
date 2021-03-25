function [] = fid_boundary(PA, EN, R2D, Z2D, FID, id, IE, CM)
% Function used to plot the FID at a particular time and position (R,Z)
% together with the particle orbit topology boundaries
%
% NB    The orbit boundaries are added manually in this scritp!
%
% INPUT
% R2D		array		radial position at which the FID has been calculated in cm
% Z2D		array		vertical position at which the FID has been calculated in cm
% PA		array		pitch angle at which the FID has been calculated 
% EN		array		energy in eV at which the FID has been calculated in cm-3 eV-1 (dW/4pi)-1
% FID		array		fast ion distribution function
% id		integer		array index of FID which is closest to the requested point in the poloidal plane
% CM 		real		Z max of the FID: 0 for auto
% IE		real		maximum injection energy for plotting

% Example
% filename = '/media/marco/WDpassport/marco/Documents/MAST/TRANSP/RUNs/29880/U16/29880U16_fi_1.cdf';
% [R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U] = fastions(filename, RR = 90, RZ = 0, CM = 0, IE = 6E4, plotyn = 1, saveyn = 0);
% fid_boundary(PA, EN, R2D, Z2D, FID, id, IE, CM);
% fid_boundary(U16pre.PA, U16pre.EN, U16pre.R2D, U16pre.Z2D, U16pre.FID, U16pre.id, U16pre.IE, U16pre.CM);



% Particle orbit Topology Boundaries 

% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s and 0.263 s, R = 0.9 m Z = 0.0 m
% ----------------------------------------------------------------------
B1.x = [0.25 0.25 0.35 0.45 0.55 0.65];
B1.y = [0 7.5 12.5 32.5 57.5 62.5]*1E3;
B2.x = [0.05 0.05];
B2.y = [0 62.5]*1E3;

% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.258 s, R = 1.1 m Z = 0.0 m
% ----------------------------------------------------------------------
B1.x = [-0.45 -0.45 -0.35 -0.35];
B1.y = [0 12.5 12.5 62.5]*1E3;
B1.y = [0 12.5 12.5 62.5];
B2.x = [-0.45 -0.45 -0.35 -0.35 -0.25 -0.25 -0.15 -0.15 -0.05 -0.05];
B2.y = [0 12.5 12.5 17.5 17.5 22.5 22.5 7.5 7.5 0]*1E3;
B2.y = [0 12.5 12.5 17.5 17.5 22.5 22.5 7.5 7.5 0];
B3.x = [-0.35 -0.35 -0.25 -0.25 -0.15 -0.15 -0.05 -0.05];
B3.y = [62.5 17.5 17.5 22.5 22.5 7.5 7.5 0]*1E3;
B3.y = [62.5 17.5 17.5 22.5 22.5 7.5 7.5 0];
B4.x = [0.05 0.05];
B4.y = [0 62.5]*1E3;
B4.y = [0 62.5];
B5.x = [0.65 0.65 0.75 0.75];
B5.y = [0 22.5 22.5 62.5]*1E3;
B5.y = [0 22.5 22.5 62.5];


% ----------------------------------------------------------------------
% Boundaries for pulse 29880 at 0.263 s, R = 1.1 m Z = 0.0 m
% ----------------------------------------------------------------------
B1.x = [-0.45 -0.45 -0.35 -0.35];
B1.y = [0 7.5 7.5 62.5];
B2.x = [-0.45 -0.45 -0.35 -0.35 -0.15 -0.15 -0.05 -0.05];
B2.y = [0 7.5 7.5 17.5 17.5 7.5 7.5 0];
B3.x = [-0.35 -0.35 -0.15 -0.15 -0.05 -0.05];
B3.y = [62.5 17.5 17.5 7.5 7.5 0];
B4.x = [0.05 0.05];
B4.y = [0 62.5];
B5.x = [0.65 0.65 0.75 0.75];
B5.y = [0 22.5 22.5 62.5];     

% Makes the plot
[pa, en] = meshgrid(PA, EN);

figure(2)
    pcolor(pa, en/1000, squeeze(FID(id,:,:))');
    shading flat
    h = xlabel('pitch angle'); set(h, 'fontsize', 12);
    h = ylabel('energy keV'); set(h, 'fontsize', 12);
    if (CM != 0)
    axis([-1 1 0 IE 0 CM])
    caxis([0 CM])
    else  
        axis([-1 1 0 IE/1E3])
    endif
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
    h = colorbar; set(h, 'title', 'cm^{-3} eV^{-1}');  
    
    % Draw the boundaries
    hold on
        h = stairs(B1.x, B1.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B2.x, B2.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B3.x, B3.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B4.x, B4.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B5.x, B5.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
    hold off
    

    
figure(3)
    pcolor(pa, en, log10(squeeze(FID(id,:,:))'));
    shading flat
    xlabel('pitch angle')
    ylabel('energy eV')
    axis([-1 1 0 6E4])
    title(['FID at R = ' num2str(R2D(id)) ' cm, Z = ' num2str(Z2D(id)) ' cm.'])
    colorbar  
    
    % Draw the boundaries
    hold on
        h = stairs(B1.x, B1.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B2.x, B2.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B3.x, B3.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B4.x, B4.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
        h = stairs(B5.x, B5.y);
        set(h, 'color', 'w', 'linestyle', '-', 'linewidth', 3)
    hold off
