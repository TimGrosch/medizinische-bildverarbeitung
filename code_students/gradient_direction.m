function [ Is ] = gradient_direction( Img )

% DESCRIPTION:
% return the absolute direction from -pi/2 to pi/2 sobel gradient
% in every point of the gradient

% convert image to double to accumolate for negative values
Img = double(Img);

% First derivative sobel mask
Dy = imfilter(Img,[1; -1],'same');
Dx = imfilter(Img,[1  -1],'same');

%Note: atan(Dy/Dx) can reach infinity if dx is zero
% Use this expression (slower but safer): 
Is = mod(atan2(Dy,Dx)+pi(), pi());

end