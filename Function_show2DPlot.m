function show2DPlot(chosen_case, chosen_layer, cases_cell)
    figure;
    imshow(squeeze(cases_cell{chosen_case}(:,:,chosen_layer)), []);
    text = append('Case: ', chosen_case, '; Layer: ', chosen_layer)
    title(text)
end