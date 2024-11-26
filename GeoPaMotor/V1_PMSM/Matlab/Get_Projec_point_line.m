function [Pp,d_min] = Projec_point_line(pl1,pl2,point)
% Coordinates of points in Input and Output are in complex

% pl1 and pl2 are two points in complex belonging to the line
% point P is the point to be projected onto the line
% Pp is the point of the projection onto the line

% Vectors P-P1 and P2-P1
Vect_Pl1_P = point - pl1;
Vect_Pl1_Pl2 = pl2 - pl1;


% Vector belonging to the line with length equal to Pp-P1
Vect_Pl1_Pp = Vect_Pl1_Pl2 * (real(Vect_Pl1_P)*real(Vect_Pl1_Pl2) + imag(Vect_Pl1_P)*imag(Vect_Pl1_Pl2)) / (abs(Vect_Pl1_Pl2)^2);

% Point Pp
Pp = pl1 + Vect_Pl1_Pp;

% distance of projection
d_min = abs(Pp-point);

end