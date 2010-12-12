function [best_cost,best_locations] = kfan_matching(costmaps, sigma, mu)

[h,w,num_k] = size(costmaps);
dist_part = 1;

dt_images = costmaps;%Something something sigma

best_cost = -1e100;
for i=1:h
    for j=1:w
        
        cost = costmaps(i,j,dist_part);
        for k = 2:num_k
            i2 = round(i+mu(k,1));
            j2 = round(j+mu(k,2));
            
            if(i2 < 1 || j2 < 1 || i2 >= h || j2>= w)
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