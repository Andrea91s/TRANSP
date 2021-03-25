% Script for the generation of the ripple UFILE

% Geoemtry
Rc = 2.09;
N = 12;

R = linspace(0.3, 1.8, 50);
Z = linspace(-1.5, 1.5, 50);
[RR, ZZ] = meshgrid(R,Z);
FF = (RR/Rc).^N;
RR1 = reshape(RR, 1, 2500);
ZZ1 = reshape(ZZ, 1, 2500);
FF1 = reshape(FF, 1, 2500);
DATA = [RR1 ZZ1 FF1];


% Open the file
filename = 'A29904.RPL';
fid = fopen (filename, 'w');


% Print header
fprintf(fid, ' 29904 MAS 2 0 6              ;-SHOT #-F(X,Y) DATA-######- 12-Nov-13\n');
fprintf(fid, '12-Nov-13                     ;-DATE - UFILE generated by MC \n');
fprintf(fid, '  0                           ;-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n');
fprintf(fid, 'R                       M     ;-INDEPENDENT VARIABLE LABEL: R-\n');
fprintf(fid, 'Z                       M     ;-INDEPENDENT VARIABLE LABEL: Z-\n');
fprintf(fid, 'RIPPLE                        ;-DEPENDENT VARIABLE LABEL-\n');                    
fprintf(fid, ' 0                            ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n');       
fprintf(fid, '      2500                    ;-# OF R PTS-\n');
fprintf(fid, '      2500                    ;-# OF Z PTS- R,Z,F(R,Z) DATA FOLLOW:\n');
fprintf(fid, '  %e %e %e %e %e %e\n', RR1); fprintf(fid, '\n');
fprintf(fid, ' %+e %+e %+e %+e %+e %+e\n', ZZ1); fprintf(fid, '\n');  
fprintf(fid, '  %e %e %e %e %e %e\n', FF1); fprintf(fid, '\n');     
fprintf(fid, ';----END-OF-DATA-----------------COMMENTS:------------------------------\n');

fclose(fid);