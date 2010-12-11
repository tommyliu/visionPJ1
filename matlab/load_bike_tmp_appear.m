b = load('../bike_tmp_appear.txt','ascii');
nparts = 6;
ndims = 16;
featsize = 50;

appearFGModel = {};

for p=1:nparts
    f = zeros(50,50,16);
    for d=1:ndims
        
        for r=1:featsize
           f(r,:,d) = b(featsize*ndims*(p-1)+featsize*(d-1)+r,:); 
        end
        
    end
    appearFGModel{p} = f;
end