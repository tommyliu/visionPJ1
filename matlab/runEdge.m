% Run edge 
%
% run_edge(imagefile, dilate, output_dir)
%   image_path    - path to image
%   output_dir    - where to save the edge map ('.' by default)
%   dilate        - dilate block size (5 by default)
%
% [NOTE]
% Code originally by Crandall (crandall@cs.cornell.edu)
% Last modification by Marynel Vazquez (marynel@cmu.edu), 12/11/10
function runEdge(image_path, output_dir, dilate)

% set default values if missing arguments
if nargin < 2, output_dir = '.'; end
if nargin < 3, dilate = 5; end

% read image
I = imread(image_path);
  
% compute edge map
Z=getSparseRepresentation(I, dilate);

% save edge map
Z=uint8(Z);
file = sprintf('%s/%s.dat', output_dir, basename(image_path));
fp=fopen(file, 'wb');
fwrite(fp, size(Z), 'integer*4');
fwrite(fp, Z, 'uint8');
fclose(fp);
  
  
