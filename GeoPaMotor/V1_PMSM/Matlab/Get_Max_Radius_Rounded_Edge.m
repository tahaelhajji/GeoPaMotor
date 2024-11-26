function Max_Radius = Get_Max_Radius_Rounded_Edge(z1,z2,z3)

% B: z1, A: z2, C:z3
AB = abs(z1-z2);
AC = abs(z2-z3);
BC = abs(z1-z3);
alpha = acos((AB^2+AC^2-BC^2)/(2*AB*AC));


Max_Radius_1 = tan(alpha/2) * AC;
Max_Radius_2 = tan(alpha/2) * AB;

Max_Radius = min([Max_Radius_1, Max_Radius_2]);


end

