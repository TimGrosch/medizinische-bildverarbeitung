% input:
% sob:      img with used sobel filter
% gpb:      img with used general probability of boundary
% can:      img with used canny filter

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


% Display the final combined edge image
figure;
title("Combined Algorithms")
imshow(final_edge,[]);

end