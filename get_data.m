% Input:
% ID:   ID of the used case

function [coronal_layer, tumor_side] = get_data(ID, data)
    
    IDs = data.CaseID;
    row = find(IDs == ID);
    column_coronal_layer = 9;
    coronal_layer = table2array(data(row, column_coronal_layer));
    column_side_tumor = 12;
    tumor_side = table2array(data(row, column_side_tumor));
    tumor_side = tumor_side{1};
end


