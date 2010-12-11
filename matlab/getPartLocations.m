% Read part locations from training file
% [partLocations,imagesPath] = getPartLocations(training_file)
%   training_file       - File with parts location per training image
%   partLocations 6x2xN - Location of each of the 6 parts in each image
%   edgesPath     {N}   - Path to edge map per image
%
% [NOTE]
% We assume the model has 6 parts.
% The training_file has the same format provided by Crandall.
function [partLocations,edgesPath] = getPartLocations(training_file)

num_parts = 6;
partLocations = [];
edgesPath = {};

imid = 1;
fp=fopen(training_file, 'r');
while(feof(fp) ~= 1)
    part_coord = zeros(6,2);
    fbasename = fgets(fp);      %basename for file
    fgets(fp);                  %file path and name
    fgets(fp);                  %likelihood
    fgets(fp);                  %just 6 2 all the time

    fprintf(fbasename);
    for i=1:num_parts
        part_coord(i,:) = str2num(fgets(fp));
    end
    fgets(fp);                  %blank line
    
    partLocations = cat(3,partLocations, part_coord);
    edgesPath{imid} = fbasename;
    
%     sprintf('%s/%s.dat', dataDirect, fbasename(1:end-2))
%     Z = readEdges(sprintf('%s/%s.dat', dataDirect, fbasename(1:end-2)));
%     
%     imagesc(Z)

    imid = imid+1;
    
end  %end of while
fclose(fp);