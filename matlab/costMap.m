% Compute the cost map for each part of the model 
% map = costMap(test_image, bgModel, fgModel)
%   test_image  MxN         - input image
%   bgModel     1xO         - background model
%   fgModel     {TxTxO}xP   - foreground model
% where O is the number of orientations considered, T is the size of the
% template and P is the number of parts. 
%
% [NOTE]
% Usually, o = 16, t = 50 and p = 6.
function the_map = costMap(test_image, bgModel, fgModel)

[M,N] = size(test_image);       % image size
num_orient = length(bgModel);   % number of orientations
num_parts = length(fgModel);    % number of parts in the model
assert(num_parts > 0);
t = size(fgModel{1},1);         % template size

% avoid 0 probabilities
bgModel = bgModel+1e-10;
for p=1:num_parts, fgModel{p} = fgModel{p}+1e-10; end

% compute logs over appereance models
bgModel = log(bgModel);

for p=1:num_parts, 
    fgModel{p}(:,:,:) = log(fgModel{p}(:,:,:));
%     fgModel{p}(:,:,1) = zeros(t,t);
end


the_map = zeros(M,N,num_parts);

tStart = tic;
for p=1:num_parts
    
    fprintf('hazah %d\n', p);
    
   for r=1:M
      for c=1:N
          
          %if (test_image(r,c) == 0), continue; end 
          
          top_row = max(1,t/2-r+1);
          bottom_row = t + min(0,(M-r)-t/2);
          left_col = max(1,t/2-c+1);
          right_col = t + min(0,(N-c)-t/2);
   
          g = 0;
          
          orientations = ...
              test_image(top_row-t/2+r:bottom_row-t/2+r, ...
                         left_col-t/2+c:right_col-t/2+c);
          
          for d = 2:16
              
              fg = fgModel{p}(top_row:bottom_row,left_col:right_col,d);
              g = g + sum(fg(orientations == d-1)) - bgModel(d)*sum(sum(orientations == d-1));
              
          end
          
          the_map(r,c,p) = g;
          
      end
   end
   
   fprintf('elapsed time: %.2f\n', toc(tStart));
   
end

% BEFORE FOR THE MAPS
% 00000000000000000000000000000000000000000000000000000000000000000000000 %
% 
% tStart = tic;
% the_map = zeros(M,N,num_parts);
% for p=1:num_parts
%     
%     fprintf('hazah %d\n', p);
%     the_map(:,:,p) = -log;
%     
%    for r=1:M
%       for c=1:N
%           
%           if (test_image(r,c) == 0), continue; end 
%           
%           top_row = max(1,t/2-r+1);
%           bottom_row = t + min(0,(M-r)-t/2);
%           left_col = max(1,t/2-c+1);
%           right_col = t + min(0,(N-c)-t/2);
%    
%           g = 0;
%           
%           orientations = ...
%               test_image(top_row-t/2+r:bottom_row-t/2+r, ...
%                          left_col-t/2+c:right_col-t/2+c);
%           
%           for d = 1:16
%               
%               fg = fgModel{p}(top_row:bottom_row,left_col:right_col,d);
%               g = g + sum(fg(orientations == d-1)) - bgModel(d)*sum(sum(orientations == d-1));
%               
%           end
%           
%           the_map(r,c,p) = g;
%           
%       end
%    end
%    
%    fprintf('elapsed time: %.2f\n', toc(tStart));
%    
% end
% 
% 00000000000000000000000000000000000000000000000000000000000000000000000 %

% % compute back probs
% fgBackProbs = {};
% for p=1:num_parts,
%    for orient=1:num_orient
%        bgOrientMinus0 = repmat(bgModel(orient)-bgModel(1),t,t);
%        assert(sum(sum(isnan(bgOrientMinus0))) == 0 && ...
%               sum(sum(isinf(bgOrientMinus0))) == 0);
%        
%        fgBackProbs{p}(:,:,orient) = ...
%            (fgModel{p}(:,:,orient) - fgModel{p}(:,:,1)) - bgOrientMinus0;
%        
%        assert(sum(sum(isnan(fgBackProbs{p}(:,:,orient)))) == 0 && ...
%               sum(sum(isinf(fgBackProbs{p}(:,:,orient)))) == 0);
%    end   
% end
% 
% the_map = zeros(M,N,num_parts);
% Z = zeros(1,num_parts);
% for p=1:num_parts
%     
%    fg = fgBackProbs{p}(:,:,1);
%    bg = repmat(bgModel(1),t,t);
%    
%    Z(p) = sum(sum(fg - bg));
%    
%    the_map(t/2+1:end-t/2,t/2+1:end-t/2,p) = repmat(Z(p),M-t,N-t);
%    assert(sum(sum(isnan(the_map(:,:,p)))) == 0);
%    
%    for r=t/2+1:M-t/2
%       for c=t/2+1:N-t/2
%          
%           this_orientation = test_image(r,c);
%           
%           % ignore no edge
%           if (this_orientation == 0), continue; end 
%           
%           top_row = max(-t/2,t/2-r)+1;
%           bottom_row = min(t/2,M-r-t/2);
%           left_col = max(-t/2,t/2-c)+1;
%           right_col = min(t/2,N-c-t/2);
%           
%           this_log_p = fgBackProbs{p}(:,:,this_orientation);
%             
%           the_map(r-top_row:r+bottom_row, c-left_col:c+right_col,p) = ...
%               conv2(the_map(r-top_row:r+bottom_row, c-left_col:c+right_col,p),...
%                     this_log_p(t/2-top_row:t-(t/2-bottom_row),t/2-left_col:t-(t/2-right_col)),'same');
%           
% %           the_map(r-top_row:r+bottom_row, c-left_col:c+right_col,p) = ...
% %               the_map(r-top_row:r+bottom_row, c-left_col:c+right_col,p) + ...
% %                     this_log_p(t/2-top_row:t-(t/2-bottom_row),t/2-left_col:t-(t/2-right_col));
% %                 
%           
%           
%       end
%    end
%    
% end



