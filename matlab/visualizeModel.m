% Visualize parts of model
% visualizePartsOfModel(image,parts,t)
%   image       - path to image
%   parts   Px2 - parts center (row,col)
%   t           - feature window size
function visualizeModel(image,parts,t)

% im = imread(image);
im = readEdges(image);

figure(); hold on;
imagesc(im);
colormap gray;
axis equal;

for p=1:size(parts,1)
    
    x = parts(p,2)-t/2;
    y = parts(p,1)-t/2;
    str = sprintf('\\fontsize{24}\\color{red}P%d',p);
    text(parts(p,2),parts(p,1),str);

    rectangle('Position', [x y t t], ...
        'EdgeColor', 'red');
end