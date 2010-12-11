% Get sparse representation 
%
% get_sparse_representation(I, edge_thresh, new_stuff, dilate_pixels)
%   I               MxN - image
%   dilate_pixels       - size of dilate block
%
% [NOTE]
% Code originally by Crandall (crandall@cs.cornell.edu)
% Last modification by Marynel Vazquez (marynel@cmu.edu), 12/11/10
function result=getSparseRepresentation(I, dilate_pixels)

% compute vertical and horizontal gradients
H=filter2([-1 0 1], I, 'same');
V=filter2([-1;0;1], I, 'same');

% compute angle of gradient
warning off MATLAB:divideByZero;
D=atan(H./V);
warning on MATLAB:divideByZero;
D(find(V==0)) = pi/2;

% groups angles by direction (1,2,3 or 4)
%   direction 1: vertical orientation
%   direction 2: south-east or north-west orientation
%   direction 3: horizontal oriantation
%   direction 4: north-east or south-west orientation
[aa bb] = size(H);
Q = zeros(aa,bb);

Q(find((D < -3*pi/8))) = 1;
Q(find((D > 3*pi/8))) = 1;
Q(find((D > -3*pi/8) .* (D < -pi/8))) = 2;
Q(find((D > -pi/8) .* (D < pi/8))) = 3;
Q(find((D > pi/8) .* (D < 3*pi/8))) = 4;

% find edges in the image
E = edge(I,'canny');

% label edges by orientation
result = Q .* E;

% distribute result in 4 different maps and dilate each map independently
result2 = zeros(size(result));
for i=1:4
result2(:,:,i) = (result == i);
result2(:,:,i) = imdilate(result2(:,:,i), ones(dilate_pixels));
end

% merge maps and label each pixel with a value in (0,15) - 2^4 posibilities
% for example, 
%   a value of 0 means no edge at any of the considered orientations
%   a value of 1 means edge matched orientation 1
result = result2(:,:,1) + result2(:,:,2)*2 + result2(:,:,3)*4 + ...
         result2(:,:,4)*8;



