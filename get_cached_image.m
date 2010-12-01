function [I]=get_cached_image(I, filter_bank)
 
[rows cols] = size(I);

% assume it's an image number

global CACHE_FILE_PREFIX 

file = sprintf('%s/images/%d.mat', CACHE_FILE_PREFIX, I);


clear result I;
load(file);

I=double(I);

