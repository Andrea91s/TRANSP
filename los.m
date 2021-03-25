function [alpha1, alpha2] = los_angles(P1, P2, P3, P4)



% Generate the line of sight
LOS = createLine (P1, P2);

% Generate the boundary circle
BC = createCircle (P3, P4);

% Calculates the intersection between the LOS
% and the boundary circle
bci = intersectLineCircle (LOS, BC);
xi1 = bci(1,1);
yi1 = bci(1,2);
xi2 = bci(2,1);
yi2 = bci(2,2);

% Calculate the angle for the two intersections
if(yi1 >= 0)
    alpha1 = atan2(yi1, xi1);
  else
    alpha1 = 2*pi+atan2(yi1, xi1);
endif

if(yi2 >= 0)
    alpha2 = atan2(yi2, xi2);
  else
    alpha2 = 2*pi+atan2(yi2, xi2);
endif


%for k =1:K
%    m1 = tan(phi(k));
%    x(p,k) = q2(p)/(m1-m2(p));
%    y(p,k) = m1*q2(p)/(m1-m2(p));
%  end