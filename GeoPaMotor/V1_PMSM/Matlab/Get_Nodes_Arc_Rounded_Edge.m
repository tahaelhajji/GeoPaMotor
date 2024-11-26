function [node1,node2,arc] = Get_Nodes_Arc_Rounded_Edge(z1,z2,z3,Radius)

% B: z1, A: z2, C:z3
AB = abs(z1-z2);
AC = abs(z2-z3);
BC = abs(z1-z3);
alpha = acos((AB^2+AC^2-BC^2)/(2*AB*AC));


node1 = z1 + (z2-z1) * (abs(z2-z1) - Radius/tan(alpha/2)) / abs(z2-z1);

node2 = z3 + (z2-z3) * (abs(z3-z2) - Radius/tan(alpha/2)) / abs(z3-z2);

arc = pi - alpha;

end

