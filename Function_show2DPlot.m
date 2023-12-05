function show2DPlot(chosen_case, chosen_layer, cases_cell, y_min, y_max)
    figure;
    imshow(squeeze(cases_cell{chosen_case}(y_min:y_max,chosen_layer,:)), []);
    %text = append('Case: ', str(chosen_case, '; Layer: ', chosen_layer)
    %title(text)
end