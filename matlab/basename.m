% Extract basename from file path
% bname = basename(filename)
%   filename - path to file
%   bname    - basename (no extension)
%
% [NOTE]
% This function is a modification of the SaliencyToolbox basename function.
% Refer to http://www.saliencytoolbox.net for more information.
function bname = basename(filename)

slash = find(filename == filesep);
if isempty(slash)
  left = 1;
else
  left = slash(end)+1;
end

thedot = find(filename == '.');
if isempty(thedot)
  right = length(filename);
else
  right = thedot(end)-1;
end

if (left > right)
  bname = filename(left:e);
else
  bname = filename(left:right);
end
