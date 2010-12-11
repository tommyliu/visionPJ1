% Read edges map
% Z = readUnknown(path)
%   path - path to edget map (.mat file generated by runEdges)
%   Z    - edges map (with values between 0-15 for each pixel)
function Z = readEdges(path)

fp = fopen(path,'r');
size = fread(fp,2,'integer*4');

Z = fread(fp,size(1)*size(2),'uint8');
Z = reshape(Z,size(1),size(2));

%imagesc(Z)