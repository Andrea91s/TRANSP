function [hi] = hintegral2d(x, y, f, a, b, c, d)
% Function that integrates the 2D histogram [x, y, f] where x and y are the bin centroids
% and f the correspoding counts in the interval [a, b] x [c, d];

% Example
% x = linspace (-0.98, 0.98, 99);   % centroids
% y = 0.2.*x.^2 + 1;                % histogram values at the centroid
% hi1 = hintegral(x, y, -5, -0.436);
% hi2 = hintegral(x, y, -0.436, 0.236);
% hi3 = hintegral(x, y, 0.236, 5);
% hi1 + hi2 + hi3
% ans =  2.1094
% sum(y)*dx
% ans =  2.1094


% Evaluates the bin width
dx = x(2) - x(1);
dy = y(2) - y(1);


% Evaluates the boundaries of the histogram bins
xb = [x(1) - dx/2 x + dx/2];
yb = [y(1) - dy/2 y + dy/2];

% 
N = length(xb);
M = length(yb);

% Find the intervals over which to integrate
[i1, i2] = find_indexes(xb, a);
if(isinf(i1) == 1)
    i1 = 1;
    i2 = 2;
    a = xb(1);
endif
[i3, i4] = find_indexes(xb, b);
if(isinf(i4) == 1)
    i3 = N-1;
    i4 = N;
    b = xb(N);
endif
[j1, j2] = find_indexes(yb, c);
if(isinf(j1) == 1)
    j1 = 1;
    j2 = 2;
    c = yb(1);
endif
[j3, j4] = find_indexes(yb, d);
if(isinf(j4) == 1)
    j3 = M-1;
    j4 = M;
    d = yb(M);
endif


% Evaluates the integral in the given interval
hi = sum(sum(f(j2:j3-1, i2:i3-1)))*dx*dy + ...
     sum(sum(f(j1, i2:i3-1)))*dx*(yb(j2)-c) + ...
     sum(sum(f(j3, i2:i3-1)))*dx*(d - yb(j3)) + ...
     sum(sum(f(j2:j3-1, i1)))*dy*(xb(i2)-a) + ...
     sum(sum(f(j2:j3-1, i3)))*dy*(b-xb(i3)) + ...
     f(j1,i1)*(xb(i2)-a)*(yb(j2)-c) + ...
     f(j1,i3)*(b-xb(i3))*(yb(j2)-c) + ...
     f(j3,i3)*(b-xb(i3))*(d-yb(j3)) + ...
     f(j3,i1)*(xb(i2)-a)*(d-yb(j3));

     %sum(f(:))*dx*dy

 
