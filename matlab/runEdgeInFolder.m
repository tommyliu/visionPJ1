% Run RunEdge for all images in specified folder
% runEdgeInFolder(folder_path, output_path)
%   images_folder     - path to images folder
%   output_path     - path to edges folder (output)
%
% [NOTE]
% A block of size 5 is used to dilate the edge maps.
function runEdgeInFolder(images_folder, output_path)

images = dir(images_folder);
for file_id=1:size(images,1)
   
    % don't process folders ni hidden files
    if (images(file_id).isdir), continue; end;
    image_path = images(file_id).name;
    if (image_path(1)=='.'), continue; end;
    display(sprintf('%s',image_path));
    
    % compute edges
    runEdge(fullfile(images_folder,image_path), output_path);

end