function [Coord_Nodes] = GetNodeCoordinate(List_Nodes)


Getelements = myfpproc.getelements;
Getvertices = myfpproc.getvertices;

Element_Nodes = Getelements(:,1:3);

Coord_Nodes = zeros(size(List_Nodes,1),2);

for k=1:length(List_Nodes)
    Node=List_Nodes(k);
    [i,j]=find(Element_Nodes==Node);
    Coord_Nodes(k,:) = Getvertices(i,[(j-1)*2+1,(j-1)*2+2]);
end


end

