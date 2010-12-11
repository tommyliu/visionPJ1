% ----------------------------------------------------------------------- %
% Configure

images_folder = '/tmpVisionPJ/test_edges/';

% ----------------------------------------------------------------------- %
% Load model data

bikeModel; 

% ----------------------------------------------------------------------- %
% Run test for given images

images = dir(images_folder);
for file_id=1:size(images,1)
   
    % don't process folders ni hidden files
    if (images(file_id).isdir), continue; end;
    image_path = images(file_id).name;
    if (image_path(1)=='.'), continue; end;
    display(image_path);
    
    % load edges image
    edges_image = readEdges(fullfile(images_folder,image_path));
    
    % detect model
    
    
    
    break;
end
