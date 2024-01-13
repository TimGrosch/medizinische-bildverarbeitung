% Input:
% img:      one layer of a volume

function[sob] = sobel_filter(img)

    sobelFilterHorizontal = [-1 0 1; -2 0 2; -1 0 1];
    sobelFilterVertical = [-1 -2 -1; 0 0 0; 1 2 1];
    
    edgesHorizontal = imfilter(double(img), sobelFilterHorizontal);
    edgesVertical = imfilter(double(img), sobelFilterVertical);
    
    sob = sqrt(edgesHorizontal.^2 + edgesVertical.^2);

end