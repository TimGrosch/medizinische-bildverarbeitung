% Input:
% V:            Volume which needs to be interpolated in z-Dimension
% data:         entry of the patient in the patienttable (one row)

function [V] = interpolate_nifti_z(V,data)

    z = data.PixelDimensionen;
    x = data.Var6;
    y = data.Var7;
    
    z_grid = 0:z:size(V,1)*z-z;
    x_grid = 0:x:size(V,2)*x-x;
    y_grid = 0:y:size(V,3)*y-y;
    
    z_query = 0:1:size(V,1)*z-z;
    x_query = x_grid;
    y_query = y_grid;
    
    [Z,X,Y] = ndgrid(z_grid,x_grid,y_grid);
    [Zq,Xq,Yq] = ndgrid(z_query,x_query,y_query);
    
    V = interpn(Z,X,Y,V,Zq,Xq,Yq);

end