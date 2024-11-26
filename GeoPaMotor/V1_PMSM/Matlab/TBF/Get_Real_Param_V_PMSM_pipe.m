function [Real_Param_pop_i] = Get_Real_Param_V_PMSM_pipe(Normal_Bd_Geo_Param_pop_i,Optim)

N_slots = Optim.Specifications.Machine.Windings.Nslots;
No_Op_Points = Optim.Specifications.Operating_Points.N_Op_Points;
Isthm = Optim.Specifications.Machine.Param.Geom.Isthm;
SO_min = Optim.Specifications.Machine.Param.Geom.SO_min;
Nphase = Optim.Specifications.Machine.Windings.Nphase;
SPP = Optim.Specifications.Machine.Windings.SPP;
p = Optim.Specifications.Machine.Arch.p;
alpha_p = pi/p;
gamma_p = pi/p/Nphase/SPP;
Fb_i_PM_min = Optim.Specifications.Machine.Param.Geom.Fb_i_PM_min;
Fb_PM_Agp_min = Optim.Specifications.Machine.Param.Geom.Fb_PM_Agp_min;
N_Geom_Param = Optim.Specifications.Machine.Param.Geom_Norm.N_Geom_Param;
Airgap_Length = Optim.Specifications.Machine.Param.Geom.Airgap_Length;

Volumic_Density_windings = Optim.Specifications.Machine.Materials.Winding.Volum_Density;

% Angle of slot pitch
Ang_Slot_Pitch = 2*pi/N_slots;

% List of Real Nodes of the rotor
Real_Nodes = zeros(1,41);
% List of Virtual Nodes of the rotor
Virtual_Nodes = zeros(1,4);

% List of Rotor Parameters
List_Rotor_Param = [];
% Real_Param_pop_i = struct([]);

for n=1:No_Op_Points
    
    Real_Param_pop_i(n).alpha_p = alpha_p;
    Real_Param_pop_i(n).gamma_p = gamma_p;
    % % Stator Parameters
    % Number of poles or magnets
    p = Optim.Specifications.Machine.Arch.p;
    
    % Speed omega in rad/s
    omega = Optim.Specifications.Operating_Points.Speed * pi / 30;
    
    % Density of PM
    rho_PM = Optim.Specifications.Machine.Materials.PM.Volum_Density;
    
    % Active Length
    Act_Leng = Normal_Bd_Geo_Param_pop_i(1)*Optim.Specifications.Machine.Param.Geom.Act_Ln_max;
    Real_Param_pop_i(n).Act_Leng = Act_Leng;
    
    % Outer Radius and Diameter of Stator
    Out_Rd_St = Normal_Bd_Geo_Param_pop_i(2)*Optim.Specifications.Machine.Param.Geom.Out_R_St_max;
    Real_Param_pop_i(n).Out_Rd_St = Out_Rd_St;
    
    % Inner Radius of Stator
    Inn_Rd_St = Normal_Bd_Geo_Param_pop_i(3)*Out_Rd_St;
    Real_Param_pop_i(n).Inn_Rd_St = Inn_Rd_St;
    
    % Height of Slot
    H_slot = Normal_Bd_Geo_Param_pop_i(4)*(Out_Rd_St-Inn_Rd_St-Isthm);
    Real_Param_pop_i(n).H_slot = H_slot;
    Real_Param_pop_i(n).Isthm = Isthm;
    
    % Slot Opening
    SO = Normal_Bd_Geo_Param_pop_i(5) * ( (pi/(p*SPP*Nphase)) * Inn_Rd_St - SO_min ) + SO_min;
    Real_Param_pop_i(n).SO = SO;
    SO_angle = SO / Inn_Rd_St;
    Real_Param_pop_i(n).SO_angle = SO_angle;
    
    % Slot Width
    SW = Normal_Bd_Geo_Param_pop_i(6) * ( (Inn_Rd_St + Isthm)*2*tan(pi/(2*p*SPP*Nphase)) - SO ) + SO;
    Real_Param_pop_i(n).SW = SW;
    SW_b = SW + 2 * tan(Ang_Slot_Pitch/2) * H_slot;
    Real_Param_pop_i(n).SW_b = SW_b;
    Real_Param_pop_i(n).Ang_Slot_Pitch = Ang_Slot_Pitch;
    
    
    
    % % Rotor Parameters
    % Outer Radius of Rotor
    Out_Rd_Rt = Inn_Rd_St - Airgap_Length;
    Real_Param_pop_i(n).Out_Rd_Rt = Out_Rd_Rt;
    
    % Param 1
    Param_1 = Normal_Bd_Geo_Param_pop_i(7) * Out_Rd_Rt;
    List_Rotor_Param = [List_Rotor_Param, Param_1];

    % Param 2
    Max_Param_2_1 = tan(alpha_p/2) * Param_1;
    Max_Param_2_2 = sqrt(Out_Rd_Rt^2-Param_1^2);
    Max_Param_2 = min([Max_Param_2_1, Max_Param_2_2]);
    Min_Param_2 = Fb_i_PM_min;
    Param_2 = Normal_Bd_Geo_Param_pop_i(8) * (Max_Param_2 - Min_Param_2) + Min_Param_2;
    List_Rotor_Param = [List_Rotor_Param, Param_2];
    
    % Param 3 : Angle of PM (half angle between magnets)
    gamma_PM = Normal_Bd_Geo_Param_pop_i(9) * pi/2;
    Param_3 = gamma_PM;
    List_Rotor_Param = [List_Rotor_Param, Param_3];
    Real_Param_pop_i(n).gamma_PM = gamma_PM;
    
    % Param 4 :
    Vnode3 = Param_1 - 1i * Param_2;
    Virtual_Nodes(3) = Vnode3;
    pt2 = Vnode3 + (cos(gamma_PM) - 1i * sin(gamma_PM));
    OC = 0 + 1i * 0;
    ptlc = Get_Intersec_Line_Cercle(Vnode3,pt2,OC,Out_Rd_Rt);
    ptc = Out_Rd_Rt * cos(alpha_p/2) - 1i * Out_Rd_Rt * sin(alpha_p/2);
    ptll = Get_Intersec_Line_Line(Vnode3,pt2,OC,ptc);
    Max_Param_4_1 = abs(ptlc - Vnode3);
    % Scalar product of vectors Vnode3-pt2 and Vnode3-ptll
    Scalar_Product = real(pt2-Vnode3) * real(ptll-Vnode3) + imag(pt2-Vnode3)*imag(ptll-Vnode3);
    if Scalar_Product > 0
        Max_Param_4_2 = abs(ptll - Vnode3);
    else
        Max_Param_4_2 = Max_Param_4_1;
    end
    Max_Param_4 = min([Max_Param_4_1, Max_Param_4_2]) - Fb_PM_Agp_min;
    Param_4 = Normal_Bd_Geo_Param_pop_i(10) * Max_Param_4;
    List_Rotor_Param = [List_Rotor_Param, Param_4];
    
    Vnode2 = Vnode3 + Param_4 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    Virtual_Nodes(2) = Vnode2;
    Real_Nodes(7) = ptc;
    Real_Nodes(5) = Out_Rd_Rt + 1i * 0;
    
    
    % Param 5
    Z1 = Vnode2 + ( cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2) );
    Z2 = Get_Intersec_Line_Line(Vnode2,Z1,OC,ptc);
    Z3 = Z2 + (-cos(gamma_PM) + 1i * sin(gamma_PM));
    ptaxe = Out_Rd_Rt + 1i *0;
    Z4 = Get_Intersec_Line_Line(Z2,Z3,OC,ptaxe);
    Min_Param_5_1 = Z4; % equal also to real(Z4)
    Min_Param_5_2 = 0;
    Min_Param_5 = max([Min_Param_5_1, Min_Param_5_2]);
    Max_Param_5 = Param_1;
    Param_5 = Normal_Bd_Geo_Param_pop_i(11) * (Max_Param_5 - Min_Param_5) + Min_Param_5;
    List_Rotor_Param = [List_Rotor_Param, Param_5];
    
    
    % Param 6
    pt1 = Param_5 + 1i *0;
    pt2 = pt1 - 1i;
    ptlc = Get_Intersec_Line_Cercle(pt1,pt2,OC,Out_Rd_Rt);
    ptll = Get_Intersec_Line_Line(pt1,pt2,OC,ptc);
    Max_Param_6_1 = abs(ptlc - pt1);
    Max_Param_6_2 = abs(ptll - pt1);
    Max_Param_6_0 = min([Max_Param_6_1, Max_Param_6_2]); % - Fb_PM_Agp_min;
    Z5 = pt1 + (0 - 1i);
    Z6 = Get_Intersec_Line_Line(pt1,Z5,Z2,Z4);
    Z7 = Get_Intersec_Line_Line(pt1,Z5,Vnode2,Z2);
    Max_Param_6_3 = abs(Z6-pt1);
    Max_Param_6_4 = abs(Z7-pt1);
    Max_Param_6 = min([Max_Param_6_0, Max_Param_6_3, Max_Param_6_4]);
    Min_Param_6_1 = 0;
    Min_Param_6_2 = Fb_i_PM_min;
    pt3 = Vnode3 + (cos(gamma_PM) - 1i * sin(gamma_PM));
    Z8 = Get_Intersec_Line_Line(pt1,pt2,Vnode3,pt3);
    Min_Param_6_3 = imag(-Z8);
    Min_Param_6 = max([Min_Param_6_1, Min_Param_6_2, Min_Param_6_3]);
    Param_6 = Normal_Bd_Geo_Param_pop_i(12) * (Max_Param_6 - Min_Param_6) + Min_Param_6;
    List_Rotor_Param = [List_Rotor_Param, Param_6];
    
    
    % Param 7
    Vnode4 = Param_5 - 1i * Param_6;
    Virtual_Nodes(4) = Vnode4;
    pt2 = Vnode4 + (cos(gamma_PM) - 1i * sin(gamma_PM));
    ptlc = Get_Intersec_Line_Cercle(Vnode4,pt2,OC,Out_Rd_Rt);
    ptll = Get_Intersec_Line_Line(Vnode4,pt2,OC,ptc);
    Max_Param_7_1 = abs(ptlc - Vnode4);
    % Scalar product of vectors Vnode4-pt2 and Vnode4-ptll
    Scalar_Product = real(pt2-Vnode4) * real(ptll-Vnode4) + imag(pt2-Vnode4)*imag(ptll-Vnode4);
    if Scalar_Product > 0
        Max_Param_7_2 = abs(ptll - Vnode4);
    else
        Max_Param_7_2 = Max_Param_7_1;
    end
    vect_perp = (Vnode4-pt2) * exp(1i * pi/2);
    Pp = Get_Intersec_Line_Line(Vnode4,pt2,Vnode3,Vnode3+vect_perp);
    x = (real(Pp-Vnode4)*real(pt2-Vnode4) + imag(Pp-Vnode4)*imag(pt2-Vnode4)) / (abs(pt2-Vnode4));
    Max_Param_7_3 = x + Param_4;
    Max_Param_7_0 = min([Max_Param_7_1, Max_Param_7_2]) - Fb_PM_Agp_min;
    Max_Param_7 = min([Max_Param_7_0, Max_Param_7_3]);
    Min_Param_7_1 = x;
    Min_Param_7_2 = 0;
    Min_Param_7 = max([Min_Param_7_1, Min_Param_7_2]);
    Param_7 = Normal_Bd_Geo_Param_pop_i(13) * (Max_Param_7 - Min_Param_7) + Min_Param_7;
    List_Rotor_Param = [List_Rotor_Param, Param_7];
    
    
    Real_Nodes(18) = Vnode3 + (Param_7 - x) * (cos(gamma_PM) - 1i * sin(gamma_PM));
    Real_Nodes(23) = Vnode4 + Param_7 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    
    
    % Param 8
    Max_Param_8 = Max_Param_7 - Param_7;
    Param_8 = Normal_Bd_Geo_Param_pop_i(14) * (Max_Param_8);
    List_Rotor_Param = [List_Rotor_Param, Param_8];
    
    
    Real_Nodes(17) = Real_Nodes(18) + Param_8 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    
    
    % Param 9
    pt1 = Param_1 - 1i * Param_2;
    pt2 = pt1 + Param_4 * (cos(gamma_PM) - 1i * sin(gamma_PM)); % replaced Param_1 with pt1
    Y = Param_5 - 1i * Param_6 + Param_7 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    vect_perp = (pt2-pt1) * exp(1i * pi/2);
    Yproj = Get_Intersec_Line_Line(pt1,pt2,Y,Y+vect_perp);
    Y1 = Yproj + Param_8 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    Y2 = Y1 + cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2);
    ptlc = Get_Intersec_Line_Cercle(Y1,Y2,OC,Out_Rd_Rt);
    ptll = Get_Intersec_Line_Line(Y1,Y2,OC,ptc);
    Max_Param_9_1 = abs(ptlc - Y1);
    Max_Param_9_2 = abs(ptll - Y1);
    ptll2 = Get_Intersec_Line_Line(Y,Yproj,OC,ptc);
    Max_Param_9_3 = abs(ptll2 - Yproj);
    Max_Param_9 = min([Max_Param_9_1, Max_Param_9_2, Max_Param_9_3]) - Fb_PM_Agp_min;
    Param_9 = Normal_Bd_Geo_Param_pop_i(15) * Max_Param_9;
    List_Rotor_Param = [List_Rotor_Param, Param_9];
    
    
    Real_Nodes(12) = Real_Nodes(17) + Param_9 * (cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2));
    
    
    
    % Param 10
    Y3 = Y1 + Param_9 * (cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2));
    Y4 = Y3 + (cos(gamma_PM) - 1i * sin(gamma_PM));
    ptlc = Get_Intersec_Line_Cercle(Y3,Y4,OC,Out_Rd_Rt);
    ptll = Get_Intersec_Line_Line(Y3,Y4,OC,ptc);
    Max_Param_10_1 = abs(ptlc - Y3);
    % Scalar product of vectors Y3-Y4 and Y3-ptll
    Scalar_Product = real(Y4-Y3) * real(ptll-Y3) + imag(Y4-Y3)*imag(ptll-Y3);
    if Scalar_Product > 0
        Max_Param_10_2 = abs(ptll - Y3);
    else
        Max_Param_10_2 = Max_Param_10_1;
    end
    Max_Param_10 = min([Max_Param_10_1, Max_Param_10_2]) - Fb_PM_Agp_min;
    Param_10 = Normal_Bd_Geo_Param_pop_i(16) * Max_Param_10;
    List_Rotor_Param = [List_Rotor_Param, Param_10];
    
    
    Vnode1 = Real_Nodes(12) + Param_10 * (cos(gamma_PM) - 1i * sin(gamma_PM));
    Virtual_Nodes(1) = Vnode1;
    
    
    
    % Param 11
    ptlc = Get_Intersec_Line_Cercle(Yproj,Y,OC,Out_Rd_Rt);
    ptll = Get_Intersec_Line_Line(Yproj,Y,OC,ptc);
    Max_Param_11_1 = abs(ptlc - Yproj) - Fb_PM_Agp_min;
    Max_Param_11_2 = abs(ptll - Yproj) - Fb_PM_Agp_min;
    Max_Param_11_3 = Max_Param_9;
    Max_Param_11 = min([Max_Param_11_1, Max_Param_11_2, Max_Param_11_3]);
    Min_Param_11_1 = abs(Y-Yproj);
    Min_Param_11_2 = Param_9;
    Min_Param_11 = max([Min_Param_11_1, Min_Param_11_2]);
    Param_11 = Normal_Bd_Geo_Param_pop_i(17) * (Max_Param_11 - Min_Param_11) + Min_Param_11;
    List_Rotor_Param = [List_Rotor_Param, Param_11];
    
    
    Virtual_Nodes(5) = Real_Nodes(18) + Param_11 * (cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2));
    Virtual_Nodes(6) = Real_Nodes(17) + Param_11 * (cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Param 12
    Param_12 = Normal_Bd_Geo_Param_pop_i(18) * ((Out_Rd_Rt-Fb_PM_Agp_min) - (Param_1+Fb_PM_Agp_min)) + (Param_1+Fb_PM_Agp_min);
    List_Rotor_Param = [List_Rotor_Param, Param_12];
    
    
    % Param 13
    Z_param12 = Param_12 + 1i * 0;
    vect_prop = (Yproj - Y1) * exp(1i * pi/2);
    point = Get_Intersec_Line_Line(Yproj,Y1,Z_param12,Z_param12+vect_prop);
    d_Z_PM = abs(point-Z_param12);
    Max_Param_13_1 = d_Z_PM;
    Max_Param_13_2 = Out_Rd_Rt - Param_12;
    Max_Param_13 = min([Max_Param_13_1, Max_Param_13_2]) - Fb_PM_Agp_min;
    Param_13 = Normal_Bd_Geo_Param_pop_i(19) * Max_Param_13;
    List_Rotor_Param = [List_Rotor_Param, Param_13];
    
    
    
    Real_Nodes(3) = Param_12 - Param_13;
    Real_Nodes(4) = Param_12 + Param_13;
    
    
    
    % Param 14
    pt1 = Param_5 - 1i * Param_6;
    pt2 = Yproj + Param_11 * (cos(gamma_PM+pi/2) - 1i * sin(gamma_PM+pi/2));
    Max_Param_14 = min([abs(pt1),abs(pt2)]);
    Param_14 = Normal_Bd_Geo_Param_pop_i(20) * Max_Param_14;
    List_Rotor_Param = [List_Rotor_Param, Param_14];
    Inn_Rd_Rt = Param_14;
    Real_Param_pop_i(n).Inn_Rd_Rt = Inn_Rd_Rt;
    
    Real_Nodes(2) = Param_14;
    Real_Nodes(6) = Param_14 * (cos(alpha_p/2) - 1i * sin(alpha_p/2));
    
    
    
    
    % Parameters and nodes for Rounded Edges
    Rad_Rd_edge_Air_min_1 = 0.2;
    % Param 15
    Max_Radius_R10_V1_V2 = Get_Max_Radius_Rounded_Edge(Real_Nodes(12),Virtual_Nodes(1),Virtual_Nodes(2));
    Rad_Rd_edge_Air_min_2 = Max_Radius_R10_V1_V2;
    Rad_Rd_edge_Air_min = min([Rad_Rd_edge_Air_min_1, Rad_Rd_edge_Air_min_2]);
    Param_15 = Normal_Bd_Geo_Param_pop_i(21) * (Max_Radius_R10_V1_V2 - Rad_Rd_edge_Air_min) + Rad_Rd_edge_Air_min;
    List_Rotor_Param = [List_Rotor_Param, Param_15];
    [RNode11,RNode12,Arc_VNode1] = Get_Nodes_Arc_Rounded_Edge(Real_Nodes(12),Virtual_Nodes(1),Virtual_Nodes(2),Param_15);
    Real_Nodes(13) = RNode11;
    Real_Nodes(14) = RNode12;
    Real_Param_pop_i(n).Arc_VNode1 = Arc_VNode1;
    
    % Param 16
    Max_Radius_V1_V2_R15 = Get_Max_Radius_Rounded_Edge(Real_Nodes(14),Virtual_Nodes(2),Real_Nodes(17));
    Rad_Rd_edge_Air_min_2 = Max_Radius_V1_V2_R15;
    Rad_Rd_edge_Air_min = min([Rad_Rd_edge_Air_min_1, Rad_Rd_edge_Air_min_2]);
    Param_16 = Normal_Bd_Geo_Param_pop_i(22) * (Max_Radius_V1_V2_R15 - Rad_Rd_edge_Air_min) + Rad_Rd_edge_Air_min;
    List_Rotor_Param = [List_Rotor_Param, Param_16];
    [RNode13,RNode14,Arc_VNode2] = Get_Nodes_Arc_Rounded_Edge(Real_Nodes(14),Virtual_Nodes(2),Real_Nodes(17),Param_16);
    Real_Nodes(15) = RNode13;
    Real_Nodes(16) = RNode14;
    Real_Param_pop_i(n).Arc_VNode2 = Arc_VNode2;
    
    % Param 17
    Max_Radius_R16_V3_V4 = Get_Max_Radius_Rounded_Edge(Real_Nodes(18),Virtual_Nodes(3),Virtual_Nodes(4));
    Rad_Rd_edge_Air_min_2 = Max_Radius_R16_V3_V4;
    Rad_Rd_edge_Air_min = min([Rad_Rd_edge_Air_min_1, Rad_Rd_edge_Air_min_2]);
    Param_17 = Normal_Bd_Geo_Param_pop_i(23) * (Max_Radius_R16_V3_V4 - Rad_Rd_edge_Air_min) + Rad_Rd_edge_Air_min;
    List_Rotor_Param = [List_Rotor_Param, Param_17];
    [RNode17,RNode18,Arc_VNode3] = Get_Nodes_Arc_Rounded_Edge(Real_Nodes(18),Virtual_Nodes(3),Virtual_Nodes(4),Param_17);
    Real_Nodes(19) = RNode17;
    Real_Nodes(20) = RNode18;
    Real_Param_pop_i(n).Arc_VNode3 = Arc_VNode3;
    
    % Param 18
    Max_Radius_V3_V4_R21 = Get_Max_Radius_Rounded_Edge(Real_Nodes(20),Virtual_Nodes(4),Real_Nodes(23));
    Rad_Rd_edge_Air_min_2 = Max_Radius_V3_V4_R21;
    Rad_Rd_edge_Air_min = min([Rad_Rd_edge_Air_min_1, Rad_Rd_edge_Air_min_2]);
    Param_18 = Normal_Bd_Geo_Param_pop_i(24) * (Max_Radius_V3_V4_R21 - Rad_Rd_edge_Air_min) + Rad_Rd_edge_Air_min;
    List_Rotor_Param = [List_Rotor_Param, Param_18];
    [RNode19,RNode20,Arc_VNode4] = Get_Nodes_Arc_Rounded_Edge(Real_Nodes(20),Virtual_Nodes(4),Real_Nodes(23),Param_18);
    Real_Nodes(21) = RNode19;
    Real_Nodes(22) = RNode20;
    Real_Param_pop_i(n).Arc_VNode4 = Arc_VNode4;
    
    
    
    % Slot Area
    Slot_Area = H_slot * (SW + SW_b)/2;

    % Fill_Factors
    Fill_Factor = Optim.Specifications.Machine.Windings.Fill_Factor;
    Real_Param_pop_i(n).Fill_Factor = Fill_Factor;
    Fill_Factor_pipe = Optim.Specifications.Machine.Windings.Fill_Factor_pipe;
    Real_Param_pop_i(n).Fill_Factor_pipe = Fill_Factor_pipe;
    
    % Param 19
    H_pipe = Normal_Bd_Geo_Param_pop_i(25) * H_slot;
    Real_Param_pop_i(n).H_pipe = H_pipe;
    
    % Evaluate the width of the pipe based on the new fill factor of the
    % slot with pipe
    % Evaluate the surface of the pipe
    S_pipe = Slot_Area * (Fill_Factor - Fill_Factor_pipe);
    % Width of the pipe
    W_pipe_max_1 = S_pipe / H_pipe;
    teta_slot = atan(((SW_b+SW)/2)/H_slot);
    W_pipe_max_2 = SW/2 + tan(teta_slot)*(H_slot-H_pipe)/2;
    W_pipe = min([W_pipe_max_1, W_pipe_max_2]);
    Real_Param_pop_i(n).W_pipe = W_pipe;
    
    
    
    % Real Nodes 8,9,10, and 11 --> apply a rounded edge to pm with radius of 0.5mm
    Rad_Rd_edge_PM_specif = 0.5;
    Max_Radius_PM_1 = Get_Max_Radius_Rounded_Edge(Real_Nodes(23),Virtual_Nodes(5),Virtual_Nodes(6));
    Max_Radius_PM_2 = Get_Max_Radius_Rounded_Edge(Virtual_Nodes(5),Virtual_Nodes(6),Real_Nodes(12));
    Rad_Rd_edge_PM = min([Rad_Rd_edge_PM_specif, Max_Radius_PM_1, Max_Radius_PM_2]);
    [RNode8,RNode9,Arc_VNode5] = Get_Nodes_Arc_Rounded_Edge(Real_Nodes(23),Virtual_Nodes(5),Virtual_Nodes(6),Rad_Rd_edge_PM);
    [RNode10,RNode11,Arc_VNode6] = Get_Nodes_Arc_Rounded_Edge(Virtual_Nodes(5),Virtual_Nodes(6),Real_Nodes(12),Rad_Rd_edge_PM);
    Real_Nodes(8) = RNode8;
    Real_Nodes(9) = RNode9;
    Real_Nodes(10) = RNode10;
    Real_Nodes(11) = RNode11;
    Real_Param_pop_i(n).Arc_VNode5 = Arc_VNode5;
    Real_Param_pop_i(n).Arc_VNode6 = Arc_VNode6;
    
    % Nodes 24 to 41 by symetry
    Real_Nodes(24:41) = conj(Real_Nodes(6:23));
    
    
    % Rotation of Nodes by alpha_p/2
    Real_Nodes = Real_Nodes * exp(1i * alpha_p/2);
    Virtual_Nodes = Virtual_Nodes * exp(1i * alpha_p/2);
    
    % Length and Thickness of PM
    L_2_PM = 2 * Param_8;
    Real_Param_pop_i(n).L_2_PM = L_2_PM;
    Th_PM = Param_11;
    Real_Param_pop_i(n).Th_PM = Th_PM;

    % Surface and Volume of all magnets
    S_2_PM = 2*p * Th_PM * L_2_PM / 1e6; % in m2
    Volume_PM = S_2_PM * Act_Leng / 1e3; % in m3
    Real_Param_pop_i(n).S_2_PM = S_2_PM;
    Real_Param_pop_i(n).Volume_PM = Volume_PM;

    % Mass of magnet
    Mass_PM = Volume_PM * rho_PM;
    Real_Param_pop_i(n).Mass_PM = Mass_PM;


    % Length of one coil
    L_Coil_act = Act_Leng;
    L_Coil_1EW = ((Inn_Rd_St + Isthm + H_slot/2) * alpha_p)/2 * pi;
    L_Coil = 2 * L_Coil_act + 2 * L_Coil_1EW;
    Real_Param_pop_i(n).L_Coil = L_Coil;

    

    % Volume of end windings (in L) (distributed windings)
    Rad_Mid_Slot = Inn_Rd_St + Isthm + H_slot/2;
    Arc_pole_Mid_Slot = Rad_Mid_Slot * alpha_p;
    H_EW = Arc_pole_Mid_Slot / 2;
    S_EW = pi * ( (Inn_Rd_St + Isthm + H_slot)^2 - (Inn_Rd_St + Isthm)^2 );
    Vol_EW = 2 * H_EW * S_EW / 1e9 * 1000;

    % Active Volume of the machine (in L)
    Vol_Act_Mach = Act_Leng * pi * (Out_Rd_St)^2 / 1e9 * 1000;

    % Volume of the machine (in L)
    Vol_Machine = Vol_Act_Mach + Vol_EW;
    Real_Param_pop_i(n).Vol_EW = Vol_EW;
    Real_Param_pop_i(n).Vol_Act_Mach = Vol_Act_Mach;
    Real_Param_pop_i(n).Volume_Machine = Vol_Machine;


    

    % Surface of one conductor turn: Supposing that a slot contains one turn
    S_turn = (Slot_Area - S_pipe) * Fill_Factor;
    Real_Param_pop_i(n).S_turn = S_turn;

    % Volume and Mass of windings
    Volume_Windings = p * SPP * Nphase * S_turn / 1e6 * L_Coil / 1e3; % in m3
    Mass_Windings = (Volume_Windings) * Volumic_Density_windings;
    Real_Param_pop_i(n).Mass_Windings = Mass_Windings;

    % Nturns
    Nturns = Optim.Specifications.Machine.Windings.Nturns;
    Real_Param_pop_i(n).Nturns = Nturns;

    % Jmax 1 from specifications
    Jmax1 = Optim.Specifications.Machine.Param.Supply.Jmax;
    % Jmax2
    Jmax2 = Jmax1; %(Optim.Specifications.Inverter.I_inv_max*Nturns)/(Slot_Area_i*Fill_Factor/Optim.Specifications.Machine.Windings.Nlayer);

    % Jmax
    Jmax = min([Jmax1, Jmax2]);

    % Current density in coil
    J = Normal_Bd_Geo_Param_pop_i(N_Geom_Param+1+(n-1)*2) * Jmax;
    Real_Param_pop_i(n).J = J;

    % Total Current per phase
    Imax = J * ((Slot_Area-S_pipe)*Fill_Factor/Optim.Specifications.Machine.Windings.Nlayer);
    Real_Param_pop_i(n).Imax = Imax;
    Real_Param_pop_i(n).Irated = Imax/sqrt(2);
    Real_Param_pop_i(n).Slot_Area = Slot_Area;

    % Control Angle
    phi0 = Normal_Bd_Geo_Param_pop_i(N_Geom_Param+2+(n-1)*2) * pi/2;
    Real_Param_pop_i(n).phi0 = phi0;

    
    %Number of phases
    Real_Param_pop_i(n).Nphase = Nphase;

    % Speed
    Speed = Optim.Specifications.Operating_Points.Speed;
    Real_Param_pop_i(n).Max_Speed = max(Speed);
    Real_Param_pop_i(n).Speed = Speed(n);
    Real_Param_pop_i(n).Speed_rad_s = Speed(n)*pi/30;

    % Number of pair poles
    Real_Param_pop_i(n).p = Optim.Specifications.Machine.Arch.p;

    % Number of slots
    Real_Param_pop_i(n).Nslots = Optim.Specifications.Machine.Windings.Nslots;

    % Number of magnets
    Real_Param_pop_i(n).N_PM = 2 * Optim.Specifications.Machine.Arch.p;

    % Remanent Flux Density (T), the used material
    Real_Param_pop_i(n).Br = Optim.Specifications.Machine.Materials.PM.Br;

    % Permeability of PM
    Real_Param_pop_i(n).mur_PM = Optim.Specifications.Machine.Materials.PM.mur;

    % Volumic density of pm
    Real_Param_pop_i(n).rho_PM = Optim.Specifications.Machine.Materials.PM.Volum_Density;

    % Volumic density of iron
    Real_Param_pop_i(n).rho_ir = Optim.Specifications.Machine.Materials.Rotor_Core.Volum_Density;

    % Permeability of iron
    Real_Param_pop_i(n).mu_iron = Optim.Specifications.Machine.Materials.Rotor_Core.mur;

    % % HB data
    % Real_Param_pop_i(n).HB_oriented = Optim.Specifications.Machine.Materials.Rotor_Core.HB_oriented;
    % Real_Param_pop_i(n).HB_perpendic = Optim.Specifications.Machine.Materials.Rotor_Core.HB_perpendic;
    % 
    % % Loss data
    % Real_Param_pop_i(n).loss = Optim.Specifications.Machine.Materials.Loss_coef;

    % Resistivity of copper
    Real_Param_pop_i(n).resistv_cu = Optim.Specifications.Machine.Materials.Winding.Resistivity;

    % Volumic Density of copper
    Real_Param_pop_i(n).rho_cu = Optim.Specifications.Machine.Materials.Winding.Volum_Density;

    % Temperature initialization
    Real_Param_pop_i(n).Temp_cu_rotor = Optim.Specifications.Machine.Thermal_Model.Rotor_Temp_init;
    Real_Param_pop_i(n).Temp_cu_stator = Optim.Specifications.Machine.Thermal_Model.Stator_Temp_init;
    
    % Number of parallel path in windings
    Real_Param_pop_i(n).N_Para_Path = Optim.Specifications.Machine.Windings.N_Para_Path;
    
    
    % Real and Virtual Nodes
    Real_Param_pop_i(n).Real_Nodes = Real_Nodes;
    Real_Param_pop_i(n).Virtual_Nodes = Virtual_Nodes;
    
    % List of Rotor Parameters
    Real_Param_pop_i(n).List_Rotor_Param = List_Rotor_Param;
    
end



end