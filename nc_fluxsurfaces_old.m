function [R, Z, RAXIS, ZAXIS, RR, ZZ] = nc_fluxsurfaces_old(filename);
% Function that returns the flux surfaces from the NetCDF file
% specified in the
%
% INPUT
% filename		NetCDF file name
%
% OUTPUT
% R				time x R x theta coordinates in cm
% Z				time x Z x theta coordinates in cm
% RAXIS
% ZAXIS
% RR				time x R x theta coordinates in cm
% ZZ				time x Z x theta coordinates in cm



% Read the moments
RMC00 = nc_read(filename, 'RMC00');
YMC00 = nc_read(filename, 'YMC00');
[M, N] = size(RMC00);
P = 3;
RMC = zeros(M,N,P);
RMS = zeros(M,N,P);
YMC = zeros(M,N,P);
YMS = zeros(M,N,P);
for p = 1:P
  if (p < 10)
    SRC = ['RMC0' num2str(p)];
    SRS = ['RMS0' num2str(p)];
    SYC = ['YMC0' num2str(p)];
    SYS = ['YMS0' num2str(p)];
  else
    SRC = ['RMC' num2str(p)];
    SRS = ['RMS' num2str(p)];
    SYC = ['YMC' num2str(p)];
    SYS = ['YMS' num2str(p)];   
  endif
  RMC(:,:,p) = nc_read(filename, SRC);
  RMS(:,:,p) = nc_read(filename, SRS);
  YMC(:,:,p) = nc_read(filename, SYC);
  YMS(:,:,p) = nc_read(filename, SYS);
end

RMC01 = nc_read(filename, 'RMC01');
RMC02 = nc_read(filename, 'RMC02');
RMC03 = nc_read(filename, 'RMC03');
RMC04 = nc_read(filename, 'RMC04');
RMC05 = nc_read(filename, 'RMC05');
RMC06 = nc_read(filename, 'RMC06');

RMS01 = nc_read(filename, 'RMS01');
RMS02 = nc_read(filename, 'RMS02');
RMS03 = nc_read(filename, 'RMS03');
RMS04 = nc_read(filename, 'RMS04');
RMS05 = nc_read(filename, 'RMS05');
RMS06 = nc_read(filename, 'RMS06');


YMC01 = nc_read(filename, 'YMC01');
YMC02 = nc_read(filename, 'YMC02');
YMC03 = nc_read(filename, 'YMC03');
YMC04 = nc_read(filename, 'YMC04');
YMC05 = nc_read(filename, 'YMC05');
YMC06 = nc_read(filename, 'YMC06');

YMS01 = nc_read(filename, 'YMS01');
YMS02 = nc_read(filename, 'YMS02');
YMS03 = nc_read(filename, 'YMS03');
YMS04 = nc_read(filename, 'YMS04');
YMS05 = nc_read(filename, 'YMS05');
YMS06 = nc_read(filename, 'YMS06');

% Read the magnetic axis position
RAXIS = nc_read(filename, 'RAXIS');
ZAXIS = nc_read(filename, 'YAXIS');

% Flux surface generation
theta = linspace(0,2*pi, 50);
[m, n] = size(RMC00);
n = n + 1;	% increase the number of the flux surfaces by 1 to include the magnetic axis
K = length(theta);
R = zeros(m, n, K);
Z = zeros(m, n, K);

% Define the cosine and sine moments
CM = zeros(P, K);
SM = zeros(P, K);
for p = 1:P
  CM(p,:) = cos(p*theta);
  SM(p,:) = sin(p*theta);
end

% Multiplies the moments and the cosines, sines
RR = zeros(M, N, K);
ZZ = zeros(M, N, K);
for m = 1:M
  RR(m,:,:) = squeeze(RMC(m,:,:))*CM + squeeze(RMS(m,:,:))*SM;
  ZZ(m,:,:) = squeeze(YMC(m,:,:))*CM + squeeze(YMS(m,:,:))*SM;
end
for n = 1:N
  for k = 1:K
    RR(:,n,k) = RMC00(:,n) + squeeze(RR(:,n,k)); 
    ZZ(:,n,k) = YMC00(:,n) + squeeze(ZZ(:,n,k)); 
  end
end

% the first flux surface contains the magnetic axis
for i = 1:m
  R(i,1,:) = RAXIS(i);
  Z(i,1,:) = ZAXIS(i);
end



%for i = 2: length(RMC00(:,1))
for i = 1: m
	%for j = 1:length(RMC00(1,:))
	for j = 2:n
		R(i,j,:) = RMC00(i,j-1) + RMC01(i,j-1)*cos(theta)   + RMS01(i,j-1)*sin(theta)   +...
				                RMC02(i,j-1)*cos(2*theta) + RMS02(i,j-1)*sin(2*theta) +...
				                RMC03(i,j-1)*cos(3*theta) + RMS03(i,j-1)*sin(3*theta) +...
				                RMC04(i,j-1)*cos(4*theta) + RMS04(i,j-1)*sin(4*theta) +...
				                RMC05(i,j-1)*cos(5*theta) + RMS05(i,j-1)*sin(5*theta) +...
						 RMC06(i,j-1)*cos(6*theta) + RMS06(i,j-1)*sin(6*theta);
		Z(i,j,:) = YMC00(i,j-1) + YMC01(i,j-1)*cos(theta)   + YMS01(i,j-1)*sin(theta)   +...
				                YMC02(i,j-1)*cos(2*theta) + YMS02(i,j-1)*sin(2*theta) +...
				                YMC03(i,j-1)*cos(3*theta) + YMS03(i,j-1)*sin(3*theta) +...
						 YMC04(i,j-1)*cos(4*theta) + YMS04(i,j-1)*sin(4*theta) +...
						 YMC05(i,j-1)*cos(5*theta) + YMS05(i,j-1)*sin(5*theta) +...
						 YMC06(i,j-1)*cos(6*theta) + YMS06(i,j-1)*sin(6*theta);
	end
end





TIME = nc_read(filename, 'TIME');
[rrt, zzt, plotindex, rts, zts, tts, rraxis, zzaxis] = fluxsurface(TIME, RR, ZZ, RAXIS, ZAXIS, 0.263, 0, 5);
[rt, zt, plotindex, rts, zts, tts, raxis, zaxis] = fluxsurface(TIME, R, Z, RAXIS, ZAXIS, 0.263, 0, 5);





figure(5, 'position', [100 100 800 800])
  hold on
    for k = 2:6:61; 
      plot(rt(k,:)'/100, zt(k,:)'/100, 'b-');
    end
    for k = 1:6:60; 
      plot(rrt(k,:)'/100, zzt(k,:)'/100, 'r-');
    end
    plot(raxis/100, zaxis/100, 'b+', rraxis/100, zzaxis/100, 'r+')
    %lev = linspace(-0.04, max(max(fs)),10);
    %contour(RR,ZZ,fs, lev, 'g');
    %plot(ras, zas, 'g+')
    %plot(RNA/100, ZNA/100, 'k+')
    xlabel('R (m)')
    ylabel('Z (m)')
    axis([0.2 1.5 -1 1])
    axis equal
   hold off
    


