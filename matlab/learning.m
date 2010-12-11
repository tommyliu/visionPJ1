dataDirect = '../train_edges';
fname = '../motorbike_training';
file = sprintf('%s.dat', fname);

% fp=fopen(file, 'r');
% 
% part_locations = [];
% while(feof(fp) ~= 1)
%     part_coord = zeros(6,2);
%     fbasename = fgets(fp);%basename for file
%     fgets(fp);%file path and name
%     fgets(fp);%likelihood
%     fgets(fp);%just 6 2 all the time
%     num_parts =6;
%     fprintf(fbasename);
%     for i=1:num_parts
%         part_coord(i,:) = str2num(fgets(fp));
%     end
%     fgets(fp);%blank line
%     
%     part_locations = cat(3,part_locations, part_coord);
% %     sprintf('%s/%s.dat', dataDirect, fbasename(1:end-2))
% %     Z = readEdges(sprintf('%s/%s.dat', dataDirect, fbasename(1:end-2)));
% %     
% %     imagesc(Z)
%     
% end  %end of while
% fclose(fp);

part_locations = getPartLocations(file);
[sigma, mu] = train_kfan(part_locations);