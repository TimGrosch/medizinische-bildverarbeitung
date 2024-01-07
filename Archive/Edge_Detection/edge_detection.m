img = rescale(slice(:,1:256));

lap = lapogau(img);

img = imgaussfilt(img,2);
% img = imadjust(img);

sob = rescale(sobel(img));
gpb = rescale(gPb(img));
can = rescale(ImprovedCanny(img,'rich'));

figure
subplot(3,1,1)
title('Sobel')
imshow(sob,[])
subplot(3,1,2)
title('Canny')
imshow(can,[])
subplot(3,1,3)
title('general Probability of Boundary')
imshow(gpb,[])

super_edge = compareEdges(sob,gpb,can);
super_edge = bwareaopen(super_edge,50);

figure
imshow(super_edge,[])


function[final_edge] = compareEdges(sob,gpb,can)

    col = size(sob,1);
    rows = size(sob,2);

    final_edge = zeros([col rows]);

    msob = rms(sob,'all');
    mcan = rms(can,'all');
    mgpb = rms(gpb,'all');

    for i = 1:col
        for j = 1:rows
            
            A = sob(i,j) >= msob;
            B = can(i,j) >= mcan;
            C = gpb(i,j) >= mgpb;
            if (A && B) || (A && C) || (B && C)
                final_edge(i,j) = sob(i,j) + can(i,j) + gpb(i,j);
            end
        end
    end

    % Rank the edge strengths at each pixel location from different algorithms
% [sorted_sobel, idx_sobel] = sort(sob(:), 'descend');
% [sorted_canny, idx_canny] = sort(can(:), 'descend');
% [sorted_gpb, idx_gpb] = sort(gpb(:), 'descend');

% Display the final combined edge image
figure;
imshow(final_edge,[]);

end

function[log] = lapogau(img)

    log_filter = ...
        [0 1 1 2 2 2 1 1 0;...
         1 2 4 5 5 5 4 2 1;...
         1 4 5 3 0 3 5 4 1;...
         2 5 3 -12 -24 -12 3 5 2;...
         2 5 0 -24 -40 -24 0 5 2;...
         2 5 3 -12 -24 -12 3 5 2;...
         1 4 5 3 0 3 5 4 1;...
         1 2 4 5 5 5 4 2 1;...
         0 1 1 2 2 2 1 1 0];

    log = imfilter(double(img),log_filter);

end

function[sob] = sobel(img)

    sobelFilterHorizontal = [-1 0 1; -2 0 2; -1 0 1];
    sobelFilterVertical = [-1 -2 -1; 0 0 0; 1 2 1];
    
    edgesHorizontal = imfilter(double(img), sobelFilterHorizontal);
    edgesVertical = imfilter(double(img), sobelFilterVertical);
    
    sob = sqrt(edgesHorizontal.^2 + edgesVertical.^2);

end