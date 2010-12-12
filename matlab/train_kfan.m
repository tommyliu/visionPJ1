function [sigma_out, mu_out] = train_kfan(part_locations)
[num_parts, dummy, num_img] = size(part_locations);

sigma = zeros(2,2,6);
mu = zeros(6,2);


best_gen_var = 1e10;
for p=1:num_parts
    distances = part_locations - repmat(part_locations(p,:,:), [num_parts,1,1]);

	gen_var_tot = 0;
    
    for i=1:num_parts
        dist_i = distances(i,:,:);
        dist_i = squeeze(dist_i)';
        %dist_i = dist_i(:,[2,1]);
        sigma(:,:,i) = cov(dist_i);
        mu(i,:) = mean(dist_i);
        gen_var_tot = gen_var_tot + det(sigma(:,:,i));
    end
    
    if (gen_var_tot < best_gen_var)
        best_gen_var = gen_var_tot;
        sigma_out = sigma;
        mu_out = mu;
    end
end