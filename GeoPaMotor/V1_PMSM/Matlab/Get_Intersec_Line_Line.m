function ptll = Get_Intersec_Line_Line(z1,z2,z3,z4)
% This function supposes that points defining the line are distinct

A = (conj(z2)*z1 - z2*conj(z1)) * (z4 - z3) - (conj(z4)*z3 - z4*conj(z3)) * (z2 - z1);
B = conj(z2 - z1) * (z4 - z3) - (z2 - z1) * conj(z4 - z3);

ptll = A  / B ;

end

