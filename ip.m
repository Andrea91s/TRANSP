function [ip] = ip(P1, P2)
% calculates the impact parameter of a line of sight


% Generate the line of sight
LOS = createLine (P1, P2);

% impact parameter
ip = distancePointLine ([0 0], LOS);
