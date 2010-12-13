addpath('../kadir-brady_scaledetector/');
imagereal = imread('/tmpVisionPJ/motorbikes/0003.jpg');
image_path = '/tmpVisionPJ/test_edges/10003.dat';

% addpath('../kadir-brady_scaledetector/');
% imagereal = imread('/tmpVisionPJ/background/image_0001.jpg');
% image_path = '/tmpVisionPJ/test_edges/770001.dat';

scale_part1 = 18;

image = readEdges(image_path);
[h,w] = size(image);

[M,N,C] = size(imagereal);
imagereal = imresize(imagereal, [h w]);

s1=5; 			%Start Scale
s2=20;			%Stop Scale
AA=0;			%Anti-aliased sampling (not available with Parzen windowing).
nbins=16;		%number of bins (set to 0 for Parzen window PDF estimation)
gsigma=1;		%sigma for Parzen window (if nbins=0. Only available on 1D)
wt=0.6;                 %threshold on Saliency values
yt=0;                   %threshold on inter-scale saliency

%1D example
div=(255/(nbins-1));	%quantisation of image.
Y=CalcScaleSaliency(uint8(double(imagereal)./div),s1,s2, nbins, gsigma,AA);
fprintf('Clustering features...');
C=GreedyCluster(Y, wt,yt);
fprintf('done\n');
figure;
imshow(imagereal);
length=size(C,2);
for i=1:length
  drawcircle(C(1,i)+1,C(2,i)+1,C(3,i)*2+1,'g');
end;

fprintf('Done with feature detection\n');

load bike_model;

dist_part = 1;
%     
% mycost = costMap(image,bg_model,fg_model);
% dt_images = distance_transform(mycost,sigma);%Something something sigma


fgModel = fg_model;
bgModel = bg_model;

% avoid 0 probabilities
bgModel = bgModel+1e-10;
fgModel{dist_part} = fgModel{dist_part}+1e-10;
% compute logs over appereance models
bgModel1 = log(bgModel);
% fgModel1= log(fgModel{dist_part}(:,:,:));
fgModel1 = fgModel{dist_part};

num_orient = size(bgModel1,2);   % number of orientations


num_feat = size(C,2);
best_cost = -1e100;

g = zeros(num_feat,1);
for i=1:num_feat

    t = size(fgModel{1},1);
    if (scale_part1 > C(3,i))
        t = floor(t/scale_part1*C(3,i));
        t = t+mod(t,2);
        fgModel_resize = zeros(t,t,16);
        for j=1:16
            fgModel_resize(:,:,j) = abs(imresize(fgModel1(:,:,j),[t t]));
        end
        fgModel_resize(:,:,:) = fgModel_resize(:,:,:)./repmat(sum(fgModel_resize(:,:,:),3),[1,1,16]);
        fgModel_resize(:,:,:) = log(fgModel_resize(:,:,:));
        
        image_resized = image;
    else
        image_resized = imresize(image,scale_part1/C(3,i));
        fgModel_resize = log(fgModel1);
 
    end

    [M,N] = size(image_resized);
    r = C(1,i);
    c = C(2,i);
    

    top_row = max(1,t/2-r+1);
    bottom_row = t + min(0,(M-r)-t/2);
    left_col = max(1,t/2-c+1);
    right_col = t + min(0,(N-c)-t/2);


    orientations = ...
      image_resized(floor(top_row-t/2+r):floor(bottom_row-t/2+r), ...
                 floor(left_col-t/2+c):floor(right_col-t/2+c));

    for d = 2:16

      fg = fgModel_resize(top_row:bottom_row,left_col:right_col,d);
      g(i) = g(i) + sum(fg(orientations == d-1)) - bgModel1(d)*sum(sum(orientations == d-1));

    end
    
end

[dummy, index] = max(g);
reference = C(1:2,index)'

scale_realz = C(3,index)


t = size(fgModel{1},1);

fgModel_resize = cell(p,1);

    if (scale_part1 > scale_realz)
        t = floor(t/scale_part1*scale_realz);
        t = t+mod(t,2);
        for p = 1:6
            fgModel_resize{p} = zeros(t,t,16);
            for j=1:16
                fgModel_resize{p}(:,:,j) = abs(imresize(fg_model{p}(:,:,j),[t t]));
            end
            fgModel_resize{p}(:,:,:) = fgModel_resize{p}(:,:,:)./repmat(sum(fgModel_resize{p}(:,:,:),3),[1,1,16]);
            %fgModel_resize{p}(:,:,:) = log(fgModel_resize{p}(:,:,:));
        end
       
        
        mu = mu./scale_part1*scale_realz;
        sigma = sigma ./ scale_part1*scale_realz;%REALLY?
        
        image_resized = image;
    else
        image_resized = imresize(image,scale_part1/scale_realz);
        mu = mu.*scale_part1/scale_realz;
        sigma = sigma .*scale_part1/scale_realz;
%         for p = 1:6
%             fgModel_resize{p}(:,:,:) = log(fg_model{p}(:,:,:));
%         end
    end

display('calling costmap');
mycost = costMap(image_resized, bg_model, fgModel_resize);

% KFAN MATCHING
% [best_cost,best_locations,dt_images] = kfan_matching(mycost, sigma, mu);
costmaps = mycost;

[h,w,num_k] = size(costmaps);
dist_part = 1;

dt_images = distance_transform(costmaps,sigma);%Something something sigma


best_cost = -1e100;
for i=max(1,reference(1)-t/2):min(reference(1)+t/2,h)
    for j=max(1,reference(2)-t/2):min(reference(2)+t/2,w)
        
        cost = costmaps(i,j,dist_part);
        
        for k = 2:num_k
            i2 = round(i+mu(k,1));
            j2 = round(j+mu(k,2));
            
            if(i2 < 1 || j2 < 1 || i2 > h || j2> w)
                cost = -1e100;
            else
                cost = cost + dt_images(i2,j2,k);
                
            end
        end
        
        if (cost > best_cost)
            best_cost = cost;
            best_ref_location = [i,j];
        end
    end
end

best_locations = repmat(best_ref_location,6,1) + mu;

% END KFAN MATCHING

close all;
figure(1);
figure(2);
for i=1:6
    a = mycost(:,:,i)';
    figure(1);subplot(2,3,i);imagesc(a); caxis([-5000 2000]); colorbar;

    figure(2);subplot(2,3,i);imagesc(dt_images(:,:,i)); colorbar;
%     figure(2);subplot(2,3,i);imagesc(cp(:,:,i)); caxis([-5000 2000]); colorbar;
end


figure(3);
imagesc(image_resized);colormap(gray);
hold on;
for i=1:6
    rectangle('Position',[best_locations(i,2),best_locations(i,1),2,2],'FaceColor','red')
end





% 
% for i=1:num_feat
%     
%     cost = mycost(C(1,i),C(2,i),dist_part);
%     
%        
%         for k = 2:num_k
%             i2 = round(C(1,i)+mu(k,1));
%             j2 = round(C(2,i)+mu(k,2));
%             
%             if(i2 < 1 || j2 < 1 || i2 > h || j2> w)
%                 cost = -1e100;
%             else
%                 cost = cost + dt_images(i2,j2,k);
%                 
%             end
%         end
%         
%         if (cost > best_cost)
%             best_cost = cost;
%             best_ref_location = [C(1,i),C(2,i)];
%         end
%     
% end
% 
% best_locations = repmat(best_ref_location,6,1) + mu;
% 
% close all;
% h1 = figure(1);
% h2 = figure(2);
% for i=1:6
%     a = mycost(:,:,i)';
%     figure(1);subplot(2,3,i);imagesc(a); caxis([-5000 2000]); colorbar; axis equal;
%     figure(2);subplot(2,3,i);imagesc(dt_images(:,:,i)); colorbar; axis equal;
% end
% 
% h3 = figure(3);
% imagesc(image);colormap(gray);
% hold on;
% for i=1:6
%     rectangle('Position',[best_locations(i,2)-5,best_locations(i,1)-5,5,5],'FaceColor','red')
% end


