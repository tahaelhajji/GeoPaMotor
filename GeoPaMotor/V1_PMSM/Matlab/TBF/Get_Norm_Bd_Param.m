function [Normal_Bd_Geo_Param_pop_i,Normal_Bd_Supp_Param_pop_i]=Get_Norm_Bd_Param(Normal_Param_pop_i,Optim)


No_Geo_Param=Optim.Specifications.Machine.Param.Geom_Norm.N_Geom_Param;
No_Op_Points=Optim.Specifications.Operating_Points.N_Op_Points;



Normal_Geo_Param_pop_i=Normal_Param_pop_i(1:No_Geo_Param,1);
Normal_Supp_Param_pop_i=Normal_Param_pop_i(No_Geo_Param+1:end,1);



Norm_Bd_List=Optim.Specifications.Machine.Param.Geom_Norm.List;
Normal_Bd_Geo_Param_pop_i=Norm_Bd_List(:,1)+(Norm_Bd_List(:,2)-Norm_Bd_List(:,1)).*Normal_Geo_Param_pop_i(:,1);



Norm_Bd_List=Optim.Specifications.Machine.Param.Supply_Norm.List;
Norm_Bd_List=repmat(Norm_Bd_List,[No_Op_Points,1]);
Normal_Bd_Supp_Param_pop_i=Norm_Bd_List(:,1)+(Norm_Bd_List(:,2)-Norm_Bd_List(:,1)).*Normal_Supp_Param_pop_i(:,1);



end