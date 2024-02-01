function[img] = connect_edges(img, threshold)
    img = double(img);

    if nargin == 1
        threshold = floor(0.075*sum(size(img))/2);
    end
    
    [row, col] = find(bwmorph(img,"endpoints"));
    conn = zeros(size(img));

    for i = 1:length(row)
        for j = 1:length(col)
            if pdist([row(i), col(i); row(j), col(j)]) <= threshold
                conn = insertShape(conn, 'Line', [col(i), row(i), col(j), row(j)], 'Color', 'white');
            end
        end
    end
    
    img = rgb2gray(rescale(img+conn));

end