% ----------------------------------------------------------------------- %
% SCRIPT CONFIGURATION 
%
% [NOTE]
% The edge maps should have been generated previously. The function runEdge 
% generates the map for a given image. The script runEdgeInFolder can be
% used to generate the edges for a bunch of images inside a folder.
% ----------------------------------------------------------------------- %

training_dat = '../motorbike_training.dat';
training_edges = '/tmpVisionPJ/train_edges';
background_edges = '/tmpVisionPJ/background_edges';
template_size = 50;

% ----------------------------------------------------------------------- %
% SPATIAL MODEL
% ----------------------------------------------------------------------- %

display('Train spatial model');
tStartS = tic;
part_locations = getPartLocations(training_dat);
[sigma, mu] = train_kfan(part_locations);
display(sprintf('(elapsed_time: %.2f)', toc(tStartS)));

% ----------------------------------------------------------------------- %
% APPEREANCE MODEL
% ----------------------------------------------------------------------- %

display('Train appereance model');
tStartA = tic;
bg_model = backgroundModel(background_edges);
fg_model = foregroundModel(training_dat, training_edges, template_size);
display(sprintf('(elapsed_time: %.2f)', toc(tStartA)));

display(sprintf('(total learning time: %.2f)', toc(tStartS)));

