function ptlc = Get_Intersec_Line_Cercle(pt1,pt2,OC,Radius)
% This function supposes that points are aligned in the following order:
% pt1 pt2 pt3 pt4 and supposes also that pt1 and pt2 are inside the cercle (OC, Radius)

vect_perp = (pt2-pt1) * exp(1i*pi/2);

pt3 = Get_Intersec_Line_Line(pt1,pt2,OC,OC+vect_perp);

dist = abs(pt3-OC);

x = sqrt(Radius^2 - dist^2);

pt4 = pt3 + x *(pt2 - pt1) / abs(pt2 - pt1);

ptlc = pt4;

end