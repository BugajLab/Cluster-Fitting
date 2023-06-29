function [point_locations] = sc_centers(image,cell_inds,cluster_size,parameter)

%Image- Image to find clusters in. %Cluster size is anticipated pixel
%diamter of a normal cluster in the image, parameter is a scalar
%multiplied by the median of the laplacian of the tophat filtered image

%Image set up
image=double(image);
sz=size(image);

%throwing away cell inds which have value of 0
cell_inds(find(image(cell_inds)==0)) = [];

%Making single cell image
min_cell = min(image(cell_inds));
cell_image = image-min_cell;
cell_image(cell_image<0) = 0;
med_cell = median(cell_image(cell_inds));
cell_image = cell_image/med_cell;

%Propagate edge points out
edge1 = getedge(cell_inds,sz,1);

%loop to creat propagated image
for i = 1:10
    
    if i == 1
    edge2 = growpoints(edge1,sz,1);
    edge2 = setdiff(edge2,cell_inds);
    else
    edge2 = growpoints(edge2,sz,1);
    edge2 = setdiff(edge2,cell_inds);
    end

%for every point in edge 2, find closest point in edge 1
e1_subs = zeros(length(edge1),2);
e2_subs = zeros(length(edge2),2);
[e1_subs(:,1),e1_subs(:,2)] = ind2sub(sz,edge1);
[e2_subs(:,1),e2_subs(:,2)] = ind2sub(sz,edge2);

e1_x = repmat(e1_subs(:,2),[1,length(e2_subs(:,2))]);
e2_x = repmat(e2_subs(:,2)',[length(e1_subs(:,2)),1]);
e1_y = repmat(e1_subs(:,1),[1,length(e2_subs(:,1))]);
e2_y = repmat(e2_subs(:,1)',[length(e1_subs(:,1)),1]);

dx = e1_x-e2_x;
dy = e1_y-e2_y;
dis_mat = (dx.^2+dy.^2).^(1/2);
[~,min_inds] = min(dis_mat);
e1_vals = cell_image(edge1);
min_vals = e1_vals(min_inds);

cell_image_save = cell_image;
cell_image = zeros(sz);
cell_image(cell_inds) = cell_image_save(cell_inds);
cell_image(edge2) = min_vals;
end


%Filter
cell_image = imgaussfilt(cell_image,0.8);

%Top hat
se = strel('disk',round(cluster_size/2));
hat = imtophat(cell_image,se);

%Gentle filter to remove high frequency noise
hat_filt = imgaussfilt(hat,0.6);

%Laplacian filter
lap_mat = -1*ones(5);
lap_mat(3,3) = 24;
lap = conv2(hat_filt,lap_mat);
lap = lap(2:sz(1)+1,2:sz(2)+1);
lap(lap<0) = 0;

%Selecting candidates
mincutoff = (parameter/5);
imspots = lap;
imspots(imspots<(mincutoff))=0;
clist=find(imspots);
clist0 = clist;
clist = intersect(clist,cell_inds);

%Removing points from same spot
    [rz,cz] = ind2sub(sz,clist);
    rzmat = repmat(rz,[1,length(clist)]); %Row subscripts of point locations
    czmat = repmat(cz,[1,length(clist)]); %Column subscripts of point locations
    rzmat2 = rzmat';
    czmat2 = czmat';
    dismat = ((rzmat-rzmat2).^2 + (czmat-czmat2).^2).^(1/2); %Distances between all points
    dismat = tril(dismat,-1); %take bottom triangle, these are the only values that matter
    dismat(find(dismat == 0)) = nan; %Makes zero values nan so they don't show up as being less than the distance cutoff
    close = find(dismat<(cluster_size-1));
    [closerz,closecz] = ind2sub(size(dismat),close); %here closerz(n) and closecz(n) are the pair of point that are too close
    p1vals = image(clist(closerz));
    p2vals = image(clist(closecz));
    valdiffs = p1vals - p2vals; %for positive values keep p1, throw out p2, and vice versa for negative values
    clist(closecz(find(valdiffs >= 0))) = nan; %keep p1, the 'row' point
    clist(closerz(find(valdiffs < 0))) = nan; %keep p2, the 'column' point
    clist(isnan(clist)) = [];

point_locations = clist;
end