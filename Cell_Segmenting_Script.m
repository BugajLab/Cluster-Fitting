%Cell Segmenting Script

%% Levy

a = 1;
x = 1:100;
b = 1;

while a == 1
    x = 1:100;
    y = (b)*x + 40 * sin(x/(10*b));
    z = cumtrapz(x,y);
    b = x(50);
end

%% Load Files

%Load files
[stack, all_filenames] = tifs4stack(); %Attached function
%[stack, all_filenames] = tifs4stack2(33); %USE for multipage tif
stack = double(stack);
sz3d = size(stack);
im_start = 1;
sz = sz3d(1:2);

if length(sz3d) == 2
    im_end = 1;
else
    im_end = sz3d(3);
end

%Variables
cell_inds_master = cell(im_end,1); %Each entry represents an image
cell_inds = cell(0); %Each entry represents a fitted cell
exclude = [];
cell_number = 0;
im_prev = [];
export_filename = [];
disp('files loaded')
%% Pick image, set lookup table
image_number = 1; %Currently displayed/analysed image
im = stack(:,:,image_number); 
im_vals = im(:);
sat_p = 1;
imin = 1 * prctile(im_vals,sat_p); %Multipliers set at 1 should give decent autoLUT, change as you want
imax = 3.5 * prctile(im_vals,100-sat_p);


%Save previous image's fits to master (if applicable)
if numel(im_prev)>0
cell_inds_master{im_prev} = cell_inds;
end
im_prev = image_number;

%Load fits from current image (if applicable)
cell_inds = cell_inds_master{image_number};
cell_number = numel(cell_inds);

%Display
    figure
    imshow(im,[imin, imax])
    if im_end == 1
        title(all_filenames)
    else
    title(all_filenames{image_number})
    end
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    for i = 1:cell_number
    [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{i});
    xc = round(mean(cell_indsx));
    yc = round(mean(cell_indsy));
    text(xc,yc,num2str(i),'Color','m','FontSize',30,'HorizontalAlignment',"center")
    end
    end
%% AUTOFIT a new cell
seg_parameter = 0.3;
if cell_number > 0
exclude = vertcat(cell_inds{:});
else
exclude = [];
end
figure('Name', 'Pick a Cell')
imshow(im,[imin,imax])
hold on
[bimage,seed] = click2seg(im,seg_parameter,exclude);
close 'Pick a Cell'

new_inds = find(bimage);
[ciy,cix] = ind2sub(sz,new_inds);                  
ymax = max(ciy);
ymin = min(ciy);
xmax = max(cix);
xmin = min(cix);
range = zeros(1,2);
range(1) = xmax-xmin;
range(2) = ymax-ymin;
xc = round(.5*range(1)+xmin);
yc = round(.5*range(2)+ymin);
range = 1.05*max(range);
shift = round(range*.5);
x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
y_limits = [max(yc-shift,1),min(yc+shift,sz(2))]; 

%Display new cell fit with red outline
figure
imshow(im,[imin, imax])
xlim(x_limits)
ylim(y_limits)
hold on
B = bwboundaries(bimage);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end
    %% more strict
    seg_parameter = seg_parameter * 1.1;
    [bimage] = segment(im,seed,seg_parameter,exclude);
    new_inds = find(bimage);

    figure
    imshow(im,[imin, imax])
    xlim(x_limits)
    ylim(y_limits)
    hold on
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    end
    %% less strict
    seg_parameter = seg_parameter * 0.9;
    [bimage] = segment(im,seed,seg_parameter,exclude);
    new_inds = find(bimage);

    figure
    imshow(im,[imin, imax])
    xlim(x_limits)
    ylim(y_limits)
    hold on
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    end
    %% keep
    new_inds_cell = cell(1,1);
    new_inds_cell{1} = new_inds;
    new_inds = new_inds_cell;
    cell_inds = [cell_inds;new_inds];
    if numel(new_inds) > 0
    cell_number = cell_number + 1;
    end
    new_inds = [];
    
    figure
    imshow(im,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    end
%% Modify a cell AUTOFITTING
        figure('Name', 'Pick a Cell to Modify')
        imshow(im,[imin,imax])
        hold on
        for i = 1:cell_number
               bimage = zeros(sz);
               bimage(cell_inds{i}) = 1;
                B = bwboundaries(bimage);
                for k = 1:length(B)
                   boundary = B{k};
                   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
                end
                [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{i});  
                xc = round(mean(cell_indsx));
                yc = round(mean(cell_indsy));
                text(xc,yc,num2str(i),'Color','r','FontSize',30,'HorizontalAlignment',"center")
        end
        [x,y] = ginput(1);
        mod_seed = [x,y];
        close 'Pick a Cell to Modify'

        x = round(x);
        y = round(y);
        a = sub2ind(sz,y,x);
        b = vertcat(cell_inds{:});
      if ismember(a,b) == 1
           c = find(b == a,1);
           cell_size = zeros(cell_number,1);
           for i = 1:cell_number
               cell_size(i) = numel(cell_inds{i});
           end
           size_cumulative = cumsum(cell_size);
           diff = size_cumulative-c;
           mod_cell = find(diff>=0,1);
      end
      seg_parameter = 0.3;
      cell_inds_not_mod_cell = cell_inds;
      cell_inds_not_mod_cell{mod_cell} = [];
      exclude = vertcat(cell_inds_not_mod_cell{:});
      mod_inds = cell_inds{mod_cell};
        [ciy,cix] = ind2sub(sz,mod_inds);   
        ymax = max(ciy);
        ymin = min(ciy);
        xmax = max(cix);
        xmin = min(cix);
        range = zeros(1,2);
        range(1) = xmax-xmin;
        range(2) = ymax-ymin;
        xc = round(.5*range(1)+xmin);
        yc = round(.5*range(2)+ymin);
        range = 1.05*max(range);
        shift = round(range*.5);
        x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
        y_limits = [max(yc-shift,1),min(yc+shift,sz(2))]; 
        
        bimage = zeros(sz);
        bimage(mod_inds) = 1;
        
        %Display new cell with red outline
            figure
            imshow(im,[imin, imax])
            xlim(x_limits)
            ylim(y_limits)
            hold on
            B = bwboundaries(bimage);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
            end
    %% more strict modification
    seg_parameter = seg_parameter * 1.1;
    [bimage] = segment(im,mod_seed,seg_parameter,exclude);
    mod_inds = find(bimage);

    figure
    imshow(im,[imin, imax])
    xlim(x_limits)
    ylim(y_limits)
    hold on
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    end  
    %% less strict modification
    seg_parameter = seg_parameter * 0.9;
    [bimage] = segment(im,mod_seed,seg_parameter,exclude);
    mod_inds = find(bimage);

    figure
    imshow(im,[imin, imax])
    xlim(x_limits)
    ylim(y_limits)
    hold on
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    end
    %% keep Autofit modification
    cell_inds{mod_cell} = mod_inds;
    mod_inds = [];

    figure
    imshow(im,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    end 
%% FREEHAND draw a new cell
            figure('Name', 'Draw a cell')
            imshow(im,[imin,imax])
            free_cell = drawfreehand;
            bimage = free_cell.createMask;
            close 'Draw a cell'
            new_inds = find(bimage);
            
            %Throw away pixels already in a cell
            if cell_number > 0
            old_cell_inds = vertcat(cell_inds{:});
            inter = ismember(new_inds,old_cell_inds);
            inter_inds = find(inter);
            new_inds(inter_inds) = [];
            end
            bimage = zeros(sz);
            bimage(new_inds) = 1;
            
            %Zoom in on current cell
            [cell_indsy,cell_indsx] = ind2sub(sz,new_inds);
            ymax = max(cell_indsy);
            ymin = min(cell_indsy);
            xmax = max(cell_indsx);
            xmin = min(cell_indsx);
            range = zeros(1,2);
            range(1) = xmax-xmin;
            range(2) = ymax-ymin;
            xc = round(.5*range(1)+xmin);
            yc = round(.5*range(2)+ymin);
            range = 1.05*max(range);
            shift = round(range*.5);
            x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
            y_limits = [max(yc-shift,1),min(yc+shift,sz(2))]; 
            
            %Display new cell with red outline
            figure
            imshow(im,[imin, imax])
            xlim(x_limits)
            ylim(y_limits)
            hold on
            B = bwboundaries(bimage);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
            end
        %% Redraw
            figure('Name', 'Draw a cell')
            imshow(im,[imin,imax])
            xlim(round(x_limits.*[0.9,1.1]))
            ylim(round(y_limits.*[0.9,1.1]))
            free_cell = drawfreehand;
            bimage = free_cell.createMask;
            close 'Draw a cell'
            new_inds = find(bimage);
            
            %Throw away pixels already in a cell
            if cell_number > 0
            old_cell_inds = vertcat(cell_inds{:});
            inter = ismember(new_inds,old_cell_inds);
            inter_inds = find(inter);
            new_inds(inter_inds) = [];
            end
            bimage = zeros(sz);
            bimage(new_inds) = 1;
            
            %Zoom in on current cell
            [cell_indsy,cell_indsx] = ind2sub(sz,new_inds);
            ymax = max(cell_indsy);
            ymin = min(cell_indsy);
            xmax = max(cell_indsx);
            xmin = min(cell_indsx);
            range = zeros(1,2);
            range(1) = xmax-xmin;
            range(2) = ymax-ymin;
            xc = round(.5*range(1)+xmin);
            yc = round(.5*range(2)+ymin);
            range = 1.05*max(range);
            shift = round(range*.5);
            x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
            y_limits = [max(yc-shift,1),min(yc+shift,sz(2))]; 
            
            %Display new cell with red outline
            figure
            imshow(im,[imin, imax])
            xlim(x_limits)
            ylim(y_limits)
            hold on
            B = bwboundaries(bimage);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
            end
        %% keep
    new_inds_cell = cell(1,1);
    new_inds_cell{1} = new_inds;
    new_inds = new_inds_cell;
    cell_inds = [cell_inds;new_inds];
    if numel(new_inds) > 0
    cell_number = cell_number + 1;
    end
    new_inds = [];
    
    figure
    imshow(im,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    end
%% Modify a cell FREEHAND
figure('Name', 'Pick a Cell to Modify')
        imshow(im,[imin,imax])
        hold on
        for i = 1:cell_number
               bimage = zeros(sz);
               bimage(cell_inds{i}) = 1;
                B = bwboundaries(bimage);
                for k = 1:length(B)
                   boundary = B{k};
                   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
                end
                [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{i});  
                xc = round(mean(cell_indsx));
                yc = round(mean(cell_indsy));
                text(xc,yc,num2str(i),'Color','r','FontSize',30,'HorizontalAlignment',"center")
        end
        [x,y] = ginput(1);
        mod_seed = [x,y];
        close 'Pick a Cell to Modify'

        x = round(x);
        y = round(y);
        a = sub2ind(sz,y,x);
        b = vertcat(cell_inds{:});
      if ismember(a,b) == 1
           c = find(b == a,1);
           cell_size = zeros(cell_number,1);
           for i = 1:cell_number
               cell_size(i) = numel(cell_inds{i});
           end
           size_cumulative = cumsum(cell_size);
           diff = size_cumulative-c;
           mod_cell = find(diff>=0,1);
      end
      
            mod_inds = cell_inds{mod_cell};   
            [cell_indsy,cell_indsx] = ind2sub(sz,mod_inds);
            ymax = max(cell_indsy);
            ymin = min(cell_indsy);
            xmax = max(cell_indsx);
            xmin = min(cell_indsx);
            range = zeros(1,2);
            range(1) = xmax-xmin;
            range(2) = ymax-ymin;
            xc = round(.5*range(1)+xmin);
            yc = round(.5*range(2)+ymin);
            range = 1.05*max(range);
            shift = round(range*.5);
            x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
            y_limits = [max(yc-shift,1),min(yc+shift,sz(2))];
      
            figure('Name', 'Redraw the cell')
            imshow(im,[imin,imax])
            xlim(round(x_limits.*[0.9,1.1]))
            ylim(round(y_limits.*[0.9,1.1]))
            free_cell = drawfreehand;
            bimage = free_cell.createMask;
            close 'Redraw the cell'
           
            %Throw away pixels already in a cell
            mod_inds = find(bimage);
            cell_inds_not_mod_cell = cell_inds;
            cell_inds_not_mod_cell{mod_cell} = [];
            old_cell_inds = vertcat(cell_inds_not_mod_cell{:});
            inter = ismember(mod_inds,old_cell_inds);
            inter_inds = find(inter);
            mod_inds(inter_inds) = [];
            
            %Zoom in on current cell
            [cell_indsy,cell_indsx] = ind2sub(sz,mod_inds);
            ymax = max(cell_indsy);
            ymin = min(cell_indsy);
            xmax = max(cell_indsx);
            xmin = min(cell_indsx);
            range = zeros(1,2);
            range(1) = xmax-xmin;
            range(2) = ymax-ymin;
            xc = round(.5*range(1)+xmin);
            yc = round(.5*range(2)+ymin);
            range = 1.05*max(range);
            shift = round(range*.5);
            x_limits = [max(xc-shift,1),min(xc+shift,sz(2))];
            y_limits = [max(yc-shift,1),min(yc+shift,sz(2))]; 
            
            %Display new cell with red outline
            figure
            imshow(im,[imin, imax])
            xlim(x_limits)
            ylim(y_limits)
            hold on
            B = bwboundaries(bimage);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
            end
        %% keep freehand modification
        cell_inds{mod_cell} = mod_inds;
        mod_inds = [];

        figure
        imshow(im,[imin, imax])
        if cell_number>=1
        hold on
            bimage = zeros(sz);
            bimage(vertcat(cell_inds{:})) = 1;
        B = bwboundaries(bimage);
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
        end 
%% Delete a cell

            figure('Name', 'Pick a Cell to Delete')
            imshow(im,[imin,imax])
            hold on
            for i = 1:cell_number
                   bimage = zeros(sz);
                   bimage(cell_inds{i}) = 1;
                    B = bwboundaries(bimage);
                    for k = 1:length(B)
                       boundary = B{k};
                       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
                    end
                    [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{i});  
                    xc = round(mean(cell_indsx));
                    yc = round(mean(cell_indsy));
                    text(xc,yc,num2str(i),'Color','r','FontSize',30,'HorizontalAlignment',"center")
            end
            [x,y] = ginput(1);
            close 'Pick a Cell to Delete'
            
            x = round(x);
            y = round(y);
            a = sub2ind(sz,y,x);
            b = vertcat(cell_inds{:});
            
            
          if ismember(a,b) == 1
               c = find(b == a,1);
               cell_size = zeros(cell_number,1);
               for i = 1:cell_number
                   cell_size(i) = numel(cell_inds{i});
               end
               size_cumulative = cumsum(cell_size);
               diff = size_cumulative-c;
               del_cell = find(diff>=0,1);
               cell_inds(del_cell) = [];
               cell_number = cell_number-1;
          end
          
              figure
                imshow(im,[imin, imax])
                if cell_number>=1
                hold on
                    bimage = zeros(sz);
                    bimage(vertcat(cell_inds{:})) = 1;
                B = bwboundaries(bimage);
                for k = 1:length(B)
                   boundary = B{k};
                   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
                end
                end
%% Pass fits from image x to image y (must be saved to master) %%
from_image = 2;
to_image = 5;

%Save current
cell_inds_master{image_number} = cell_inds;

%Pass
cell_inds_master{to_image} = cell_inds_master{from_image};
cell_inds = cell_inds_master{image_number};
cell_number = numel(cell_inds);
disp('done')
%% Import fits from file
          [cell_filename,~]= uigetfile('*.mat');
          cim = load(cell_filename);
          cell_inds_master = cim.cell_inds_master;
          cell_inds = cell_inds_master{image_number};
          cell_number = numel(cell_inds);
%% Save fits to file

%Update master array
if cell_number>0
    cell_inds_master{im_prev} = cell_inds;
end

%Choose file location
export_path = uigetdir;
location_name = [export_path '\Cell_indices.mat'];            
save(location_name, 'cell_inds_master');
           
%% Display cell fits
    figure
    imshow(im,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    for i = 1:cell_number
    [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{i});
    xc = round(mean(cell_indsx));
    yc = round(mean(cell_indsy));
    text(xc,yc,num2str(i),'Color','m','FontSize',30,'HorizontalAlignment',"center")
    end
    end
    
%% Clusters and Centers on current image (just for choosing parameters)
   center_parameter_v = 35;
   cluster_parameter_v = 3.5;
   cluster_size = 5;
   
   centers_v = cell(cell_number,1);
   clusters_v = cell(cell_number,1);
   

   for i = 1:cell_number
       centers_v{i} = sc_centers(im,cell_inds{i},cluster_size,center_parameter_v);
       centers_v_i = centers_v{i};
       clusters_v_i = cell(numel(centers_v_i),1);
       for j = 1:numel(centers_v_i)
           centers_v_ij = centers_v_i(j);
                    if j == 1
                        ex = centers_v_i;
                        ex(j) = [];
                    else
                        ex1 = centers_v_i;
                        ex1(j) = [];
                        ex2 = vertcat(clusters_v_i{1:(j-1)});
                        ex = cat(1,ex1,ex2);
                    end
       [clusters_v_i{j},new_center] = segment_sc_clump2(im,cell_inds{i},centers_v_ij,cluster_parameter_v,ex);
       centers_v_i(j) = new_center;
       end
       clusters_v{i} = clusters_v_i;
       centers_v{i} = centers_v_i;
   end
 
im_i = im;
dis_cell = cell_inds;
dis_centers = centers_v;
dis_clusters = clusters_v;

    figure
    imshow(im_i,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(dis_cell{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    end
   
   
   %Cluster Centers
            center_inds = vertcat(dis_centers{:});
            [ccy,ccx] = ind2sub(sz,center_inds);
            plot(ccx,ccy,'g.')
            hold on
            
            
            %Outlines
            bimage = zeros(sz);
                    cluster_inds = vertcat(dis_clusters{:});
                    cluster_inds = vertcat(cluster_inds{:});
                    bimage(cluster_inds) = 1;
                    B = bwboundaries(bimage);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:,2), boundary(:,1), 'b.')
                    end

%% Batch Centers and Clusters
   center_parameter = 35;
   cluster_parameter = 3.5;
   
   center_master = cell(im_end,1);
   cluster_master = cell(im_end,1);
   cluster_size = 5;
   
   cell_inds_master{image_number} = cell_inds;
   
   for i = 1:im_end %image loop
       im_i = stack(:,:,i);
       cell_inds_i = cell_inds_master{i};
       cell_n = numel(cell_inds_i);
       center_master_i = cell(cell_n,1);
       cluster_master_i = cell(cell_n,1);
       
       for j = 1:cell_n %cell loop
           
           %Cluster centers
           cell_inds_ij = cell_inds_i{j};
           center_master_ij = sc_centers(im_i,cell_inds_ij,cluster_size,center_parameter);
           center_master_i{j} = center_master_ij;
           cluster_number = numel(center_master_ij);
           cluster_master_ij = cell(cluster_number,1);
           
           %Cluster Segmentation
           for k = 1:cluster_number %cluster loop
               center_master_ijk = center_master_ij(k);
                    if k == 1
                        ex = center_master_ij;
                        ex(k) = [];
                    else
                        ex1 = center_master_ij;
                        ex1(k) = [];
                        ex2 = vertcat(cluster_master_ij{1:(k-1)});
                        ex = cat(1,ex1,ex2);
                    end
                    [cluster_master_ijk,new_center] = segment_sc_clump2(im_i,cell_inds_ij,center_master_ijk,cluster_parameter,ex);
                    cluster_master_ij{k} = cluster_master_ijk;
                    center_master_ij(k) = new_center;
           end
           cluster_master_i{j} = cluster_master_ij;
           center_master_i{j} = center_master_ij;
       end
       center_master{i} = center_master_i;
       cluster_master{i} = cluster_master_i;
       
       dis1 = 'Image';
       dis2 = num2str(i);

       disp([dis1 dis2])
   end
    
%% Cluster and Outlines display

i = 3; %image which will be displayed

im_i = stack(:,:,i);
dis_cell = cell_inds_master{i};
dis_centers = center_master{i};
dis_clusters = cluster_master{i};

    figure
    imshow(im_i,[imin, imax])
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(dis_cell{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    end
   
   
   %Cluster Centers
            center_inds = vertcat(dis_centers{:});
            [ccy,ccx] = ind2sub(sz,center_inds);
            plot(ccx,ccy,'g.')
            hold on
            
            
            %Outlines
            bimage = zeros(sz);
                    cluster_inds = vertcat(dis_clusters{:});
                    cluster_inds = vertcat(cluster_inds{:});
                    bimage(cluster_inds) = 1;
                    B = bwboundaries(bimage);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:,2), boundary(:,1), 'b.')
                    end
       
%% Data Export & masks

        export_filename = 'test';
        
if isempty(export_filename) == 1
   disp('Specify a filename friend')
else
   %Cells
            image_number = []; %Table entry specifying which image each pixel came from
            cell_number = []; %Table entry specifying which cell each pixel came from
            x = []; %Table entries specifying pixel location
            y = [];
            intensity = []; %Table entry specifiying intensity of each pixel
            total_cell_number = []; %Table entry with cell number from whole image set
            running_cell_count = 0;
            
            for i = 1:im_end %looping over images
                 im_i = stack(:,:,i);
                 cell_inds_master_i =  cell_inds_master{i}; %Each entry is an image -> each entry is a cell, entire array represents one image
                 if isempty(cell_inds_master_i) == 1
                    continue
                 end
                 szcmi = size(cell_inds_master_i);
                 image_number0 = i*ones(length(vertcat(cell_inds_master_i{:})),1);
                 image_number = [image_number ; image_number0]; %Table entry specifying which image each pixel came from
                 for j = 1:szcmi(1) %looping over cells
                     running_cell_count = running_cell_count+1; %for total cell number
                     cell_inds_master_ij = cell_inds_master_i{j}; %each entry representing image -> each entry representing cell
                     if isempty(cell_inds_master_ij) == 1
                        continue
                     end
                     szcmij = size(cell_inds_master_ij);
                     cell_number0 = j*ones(szcmij(1),1); %each entry representing cell -> each entry representing pixel
                     total_cell_number0 = running_cell_count*ones(szcmij(1),1);
                     cell_number = [cell_number ; cell_number0]; %Table entry specifying which cell each pixel came from
                     total_cell_number = [total_cell_number;total_cell_number0];
                     [y0,x0] = ind2sub(sz,cell_inds_master_ij(:)); 
                     x = [x ; x0]; %Table entries specifying pixel location
                     y = [y ; y0];
                     intensity0 = im_i(cell_inds_master_ij);
                     intensity = [intensity ; intensity0]; %Table entry specifing pixel intensity
                 end
            end
            T1 = table(x,y,intensity,image_number,cell_number,total_cell_number);
            
      %Clusters
            image_number = []; %Table entry specifying which image each pixel came from
            cell_number = []; %Table entry specifying which cell each pixel came from
            total_cell_number = [];
            x = []; %Table entries specifying pixel location
            y = [];
            intensity = []; %Table entry specifiying intensity of each pixel 
            
            cluster_number = [];
            running_cell_count = 0;
               
            for i = 1:im_end %looping over images
                im_i = stack(:,:,i);
                cluster_master_i = cluster_master{i}; %Each entry is an image -> each entry is a cell
                if isempty(cluster_master_i) == 1
                   continue
                end
                reducingi = vertcat(cluster_master_i{:}); %Each entry is a cell - > each entry is a clump
                reducingi = vertcat(reducingi{:}); %Each entry is a cell - > each entry is a pixel
                image_number0 = i*ones(length(reducingi),1);
                image_number = [image_number ; image_number0]; %Table entry specifying which image each pixel came from
                szcmi = size(cluster_master_i);
                for j = 1:szcmi(1) %looping over cells
                    running_cell_count = running_cell_count+1;
                    cluster_master_ij = cluster_master_i{j}; %Each entry is a cell -> each entry is a clump
                    if isempty(cluster_master_ij) == 1
                       continue
                    end
                    reducingj = vertcat(cluster_master_ij{:}); %Each entry is a clump -> each entry is a pixel
                    szrj = size(reducingj);
                    cell_number0 = j*ones(szrj(1),1); 
                    cell_number = [cell_number ; cell_number0]; %Table entry specifying which cell each pixel came from
                    total_cell_number0 = running_cell_count*ones(szrj(1),1); 
                    total_cell_number = [total_cell_number ; total_cell_number0];
                    szcmij = size(cluster_master_ij);
                    for k = 1:szcmij(1) %looping over clumps
                        cluster_master_ijk = cluster_master_ij{k}; %Each entry is clump -> each entry is a pixel
                        if isempty(cluster_master_ijk) == 1
                           continue
                        end
                        szcmijk = size(cluster_master_ijk);
                        cluster_number0 = k*ones(szcmijk(1),1);
                        cluster_number = [cluster_number ; cluster_number0]; %Table entry specifying which clump each pixel came from
                        [y0,x0] = ind2sub(sz,cluster_master_ijk(:)); 
                        x = [x ; x0]; %Table entries specifying pixel location
                        y = [y ; y0];
                        intensity0 = im_i(cluster_master_ijk(:));
                        intensity = [intensity ; intensity0]; %Table entry specifiying intensity of each pixel
                    end
                end
            end
          
            T2 = table(x,y,intensity,image_number,cell_number,cluster_number,total_cell_number);
        
        %Information
            Image_Number = [1:im_end]';
            File_Name = all_filenames;
            sz_in = size(Image_Number);
            if sz_in(1) == 1
                T3 = cell(4,1);
                T3{1} = Image_Number;
                T3{2} = File_Name;
                T3{3} = center_parameter*ones(im_end,1);
                T3{4} = cluster_parameter*ones(im_end,1);
                T3 = cell2table(T3);
            else
             T3 = table(Image_Number,File_Name,center_parameter*ones(im_end,1),cluster_parameter*ones(im_end,1));
            end
            
        %Making masks   
        mask = zeros(sz(1),sz(2),im_end);
        
        for i = 1:im_end
            cell_inds_i = cell_inds_master{i};
            if isempty(cell_inds_i) == 0
            center_inds_i = center_master{i};
            cluster_inds_i = cluster_master{i};
            cluster_inds_i = vertcat(cluster_inds_i{:});
            mask_i = zeros(sz(1),sz(2));
            mask_i(vertcat(cell_inds_i{:})) = 1;
            if isempty(cluster_inds_i) == 0
            mask_i(vertcat(cluster_inds_i{:})) = 2;
            end
            mask_i(vertcat(center_inds_i{:})) = 3;
            mask(:,:,i) = mask_i;
            end
        end
        
        %Centers
        
        
        %Writing files
        spreadsheet_export_path = uigetdir;
        spreadsheet_name = [spreadsheet_export_path '\' export_filename];
        folder = [spreadsheet_export_path '\' export_filename];
        mkdir(folder);
        
        filename_cell = [folder '\' export_filename '_cell' '.csv'];
        filename_clusters = [folder '\' export_filename '_clusters' '.csv'];
        filename_info = [folder '\' export_filename '_info' '.csv'];
        
        writetable(T1,filename_cell)
        writetable(T2,filename_clusters)
        writetable(T3,filename_info)
        
        %Writing masks
        mask_path = [spreadsheet_export_path '\' export_filename '_Masks'];
        mkdir(mask_path);
            if isa(all_filenames,'cell') == 1
                for i = 1:im_end-im_start + 1
                    name = [mask_path '\' all_filenames{i} '_mask.tif'];
                    imwrite(uint16(mask(:,:,i)), name,'tif');
                end
            elseif isa(all_filenames,'char') == 1
                    name = [mask_path '\' all_filenames '_mask.tif'];
                    imwrite(uint16(mask(:,:,i)), name,'tif');
            end
   
end   
  disp('All Done!') 
   
%% Export PNGs showing cell fits
%Save previous image's fits to master (if applicable)
if numel(im_prev)>0
cell_inds_master{im_prev} = cell_inds;
end
im_prev = image_number;

png_export_path = uigetdir;
for i = 1:im_end
    
    %calculate lut
    im = stack(:,:,i); 
    im_vals = im(:);
    sat_p = 1;
    imin = 1 * prctile(im_vals,sat_p); %Multipliers set at 1 should give decent autoLUT, change as you want
    imax = 1 * prctile(im_vals,100-sat_p);
    
    if im_end == 1
        name = all_filenames(1);
        name = [name,'.png'];
    else
    name = all_filenames{i};
    name = [name '.png'];
    end
    
    cell_inds = cell_inds_master{i};
    cell_number = numel(cell_inds);
    fig = imshow(im, [imin,imax]);
    
    if cell_number>=1
    hold on
        bimage = zeros(sz);
        bimage(vertcat(cell_inds{:})) = 1;
    B = bwboundaries(bimage);
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    for j = 1:cell_number
    [cell_indsy,cell_indsx] = ind2sub(sz,cell_inds{j});
    xc = round(mean(cell_indsx));
    yc = round(mean(cell_indsy));
    text(xc,yc,num2str(j),'Color','m','FontSize',30,'HorizontalAlignment',"center")
    end
    end
    
    
    saveas(gcf,[png_export_path name]);
    close all
    
end


%% Background Subtract

for i=1:im_end
    
    c_image = stack(:,:,i);
    c_image = c_image-min(c_image(:));
    c_image(c_image<0) = 0;
    stack(:,:,i) = c_image;
    
    
end

disp('donzo')


   

   

   