function [gPb_response] = gPb(img)

% Calculate gradients using Sobel operators
[Gx, Gy] = gradient(img);

% Compute gradient magnitude and orientation
grad_mag = sqrt(Gx.^2 + Gy.^2);
grad_orient = atan2(Gy, Gx);

% Compute gPb response using oriented filters
gPb_response = compute_gPb(grad_mag, grad_orient); % Function to compute gPb

% Normalize gPb response to [0, 1] for visualization purposes
gPb_response = (gPb_response - min(gPb_response(:))) / (max(gPb_response(:)) - min(gPb_response(:)));

end

% ucm_map = contours2ucm(gPb_response, 'imageSize'); % Use appropriate function for hierarchical segmentation
% 
% % Normalize ucm_map for visualization (if necessary)
% ucm_map = (ucm_map - min(ucm_map(:))) / (max(ucm_map(:)) - min(ucm_map(:)));
% 
% % Display the ultrametric contour map
% figure;
% imshow(ucm_map);

function gPb_response = compute_gPb(grad_mag, grad_orient)
    % Define orientations for Gabor filters
    num_orientations = 16;  % Number of orientations
    orientations = linspace(0, pi, num_orientations + 1);
    orientations = orientations(1:end-1);
    
    % Initialize gPb response matrix
    [rows, cols] = size(grad_mag);
    gPb_response = zeros(rows, cols);
    
    % Apply Gabor filters at different orientations
    for i = 1:num_orientations
        % Create Gabor filter for each orientation
        lambda = 3;  % Gabor wavelength
        sigma = 1;   % Gabor standard deviation
        filter = gabor_filter(lambda, sigma, orientations(i));  % Function to generate Gabor filter
        
        % Convolve gradient magnitude with Gabor filter
        response = conv2(grad_mag, filter, 'same');
        
        % Accumulate responses across orientations
        gPb_response = gPb_response + response;
    end
    
    % Normalize gPb_response to [0, 1] for visualization
    gPb_response = (gPb_response - min(gPb_response(:))) / (max(gPb_response(:)) - min(gPb_response(:)));
end

% Function to generate Gabor filter
function filter = gabor_filter(lambda, sigma, orientation)
    % Define Gabor filter parameters
    k = 2 * pi / lambda;
    theta = orientation;
    
    % Create Gabor filter
    sz = 3 * ceil(sigma);
    [x, y] = meshgrid(-sz:sz);
    x_theta = x * cos(theta) + y * sin(theta);
    y_theta = -x * sin(theta) + y * cos(theta);
    filter = exp(-(x_theta.^2 + y_theta.^2) / (2 * sigma^2)) .* cos(k * x_theta);
end
