function stuff = kfan_matching(costmaps, sigma, mu)

[h,w,num_k] = size(costmaps);
dist_part = 1;


for i=1:h
    for j=1:w
        
        cost = costmaps(i,j,dist_part);
        for k = 1:num_k
            
            
        end
        
        
    end
end
