% Generate foreground appereance model
%
% model = backgroundModel(training_file, edges_folder, t)
%   training_file            - File with parts location per training image
%   edges_folder             - path to background edges folder
%   t                        - size of template
%   model       {txtx16}x6 - foreground model
%
% [NOTE]
% The model is a cell of 6 elements (6 parts). Each element is a matrix of
% txtx16 elements, where for each template we have probabilities for all
% possible edge orientations (from 0 to 15)
function model = foregroundModel(training_file, edges_folder, t)

model = {};
num_parts = 6;

% get part locations
[part_locations, edges_file] = getPartLocations(training_file);
% edges_file

% process each image



