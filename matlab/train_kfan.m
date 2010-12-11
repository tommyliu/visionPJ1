function [sigma, mu] = train_kfan(part_locations)
[num_parts, dummy, num_img] = size(part_locations);

sigma_ref = [0,0,0,0];
mu_ref = [0,0];

sigma = zeros(2,2,6);
mu = zeros(6,2);


distances = part_locations - repmat(part_locations(1,:,:), [num_parts,1,1]);

for i=1:num_parts
    dist_i = distances(i,:,:);
    dist_i = squeeze(dist_i)';
    sigma(:,:,i) = cov(dist_i);
    mu(i,:) = mean(dist_i);
end

