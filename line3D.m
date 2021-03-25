function [x, y, z, t] = line3D(P1, P2, N)
% Function that calculates the straight line in 3D that connects
% points P1 and P2
%
% INPUT
% P1		point 1 coordinates
% P2		point 2 coordinates
% N			number of points between P1 and P2
%
% OUTPUT
% x, y, z	line evaluation
% t			distance coordinate along the line joining P1 and P2
%           such that t(1) = 0 at P1 and t(end) = d at P2

% Points coordinates
x1 = P1(1); y1 = P1(2); z1 = P1(3);
x2 = P2(1); y2 = P2(2); z2 = P2(3);


% Distance between points
d = sqrt((x2-x1)^2 + (y2-y1).^2 + (z2-z1)^2);


% Direction cosines
l = (x2-x1)/d;
m = (y2-y1)/d;
n = (z2-z1)/d;


% Equation of line joining P1 and P2 in parametric form
t = linspace(0, d, N);
x = x1 + l*t;
y = y1 + m*t;
z = z1 + n*t;



