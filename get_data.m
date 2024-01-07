% Input:
% ID:   ID of the used case

function [coronal_layer] = get_data(ID, data)
    
    IDs = data.CaseID;
    row = find(IDs == ID);
    column_coronal_layer = 9;
    coronal_layer = table2array(data(row, column_coronal_layer));


end


