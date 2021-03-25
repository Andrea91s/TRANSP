% Analysis of the fast ion distribution before and after the sawtooth
% for different modelling of the sawtooth


filename = '/home/andrea/Documents/TRANSP/29880/29880U16_fi_3.cdf';
[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U16] = fastions(filename, RR = 110, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
U16.filename = filename;
clear -x U*;

filename = '/home/andrea/Documents/TRANSP/29880/29880U16_fi_3.cdf';
[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U31] = fastions(filename, RR = 110, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
U31.filename = filename;
clear -x U*;

filename = '/home/andrea/Documents/TRANSP/29880/29880U16_fi_3.cdf';
[R2D, Z2D, PA, EN, FID, BMVOL, X, Y, FIDD, FIDT, FIDDT, FIDN, BDENS2T, FIDE, FIDP, r, er, id, U32] = fastions(filename, RR = 110, RZ = 0, CM = 0, IE = 6E4, plotyn = 0, saveyn = 0);
U32.filename = filename;
clear -x U*;

% Does some plots
    
figure(2)
    h = plot(U16.r, U16.er, 'b', U31.r, U31.er, 'r', U32.r, U32.er, 'g'); 
    set (h, 'linewidth', 2);
    xlabel ('R (cm)')
    ylabel ('Emissivity (cm^{-3} s^{-1})')
    legend('U16 - Standard Kadomstev', 'U31 - Kadomstev mixing', 'U32 - Porcelli Mixing')
    title(strrep(U16.filename, '_', ' - '))
  return 
 
load DIFFmap.mat 

figure(1, 'position', [100 100 1600 600])
  subplot(1,3,1)
    pcolor(U31.X, U31.Y, (U16.FIDD - U31.FIDD)/max(max(U16.FIDD))); shading interp; colormap(DIFFmap); axis tight; axis equal; caxis([-0.5 0.5]); colorbar
  subplot(1,3,2)
    pcolor(U31.X, U31.Y, (U16.FIDD - U32.FIDD)/max(max(U16.FIDD))); shading interp; colormap(DIFFmap); axis tight; axis equal; caxis([-0.5 0.5]); colorbar
  subplot(1,3,3)
    pcolor(U31.X, U31.Y, (U31.FIDD - U32.FIDD)/max(max(U31.FIDD))); shading interp; colormap(DIFFmap); axis tight; axis equal; caxis([-0.1 0.1]); colorbar    
    xlabel ('R (cm)')
    ylabel ('Z (cm)')
    title ('(U31.FIDD - U32.FIDD)/max(max(U31.FIDD))')
    
 



