function [List_dist] = Get_Proj_List(List_Z,Z1,Z2)

% Coordinates of List_Z point in Input and Output are in complex

% pl1 and pl2 are two List_Zs in complex belonging to the line
% List_Z P is the List_Z to be projected onto the line
% Pp is the List_Z of the projection onto the line
 
% Vectors P-P1 and P2-P1
Vect_Pl1_P = List_Z - Z1;
Vect_Pl1_Pl2 = Z2 - Z1;


% Vector belonging to the line with length equal to Pp-P1
Vect_Pl1_Pp = Vect_Pl1_Pl2 .* (real(Vect_Pl1_P).*real(Vect_Pl1_Pl2) + imag(Vect_Pl1_P).*imag(Vect_Pl1_Pl2)) ./ (abs(Vect_Pl1_Pl2).^2);

% List_Z Pp
Pp = Z1 + Vect_Pl1_Pp;

% distance of projection
List_dist = abs(Pp-List_Z);

end

