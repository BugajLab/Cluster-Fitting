function [edge_indices] = getedge(indices,image_size,number_of_iterations)

sz = image_size;
n = number_of_iterations;
inds = indices;

bimage = zeros(sz);
bimage(inds) = 1;
%bimage = imfill(bimage,'holes');

b_inds = [];

for i = 1:n

b = bwboundaries(bimage);
b = vertcat(b{:});
b_inds = [b_inds ; sub2ind(sz,b(:,1),b(:,2))];
inter = ismember(inds,b_inds);
inter_inds = find(inter);
inds(inter_inds) = [];
bimage = zeros(sz);
bimage(inds) = 1;

end

edge_indices = b_inds;