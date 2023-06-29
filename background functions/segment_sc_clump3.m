function [clump_inds,new_center] = segment_sc_clump3(image,cell_inds,cluster_center,parameter,exclude)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%See if we need to shift the center a little
sz = size(image);
[y,x] = ind2sub(sz,cluster_center);
top=y-2;
bottom=y+2;
left=x-2;
right=x+2;
        if top<1 
           top=1;
        end
        if bottom>sz(1)
           bottom=sz(1);
        end
        if left<1
            left=1;
        end
        if right>sz(2)
           right=sz(2);
        end
local_area=image(top:bottom,left:right);
[~,new_max_ind] = max(local_area,[],'all','linear');
[new_max_y,new_max_x] = ind2sub(size(local_area),new_max_ind);
x = new_max_x + left - 1;
y = new_max_y + top - 1; 
temp_center = sub2ind(sz,y,x);

if isempty(intersect(temp_center,cell_inds)) == 0
    cluster_center = temp_center;
else
end

cutoff = parameter;
subs = zeros(1,2);
[subs(:,1),subs(:,2)] = ind2sub(sz,cluster_center);
bimage(cluster_center) = 1;
grow = 1;
subs_new = zeros(1,2,8);
rejects = exclude;
cutoff_high = 2*image(cluster_center);
iterations = 1;

     while grow == 1

         
        %subscripts of neighboring points
        subs_new(:,1,1) = subs(:,1)+1;
        subs_new(:,2,1) = subs(:,2);
        
        subs_new(:,1,2) = subs(:,1);
        subs_new(:,2,2) = subs(:,2)+1;
        
        subs_new(:,1,3) = subs(:,1)-1;
        subs_new(:,2,3) = subs(:,2);
        
        subs_new(:,1,4) = subs(:,1);
        subs_new(:,2,4) = subs(:,2)-1;
        
        subs_new(:,1,5) = subs(:,1)+1;
        subs_new(:,2,5) = subs(:,2)+1;
        
        subs_new(:,1,6) = subs(:,1)-1;
        subs_new(:,2,6) = subs(:,2)-1;
        
        subs_new(:,1,7) = subs(:,1)+1;
        subs_new(:,2,7) = subs(:,2)-1;
        
        subs_new(:,1,8) = subs(:,1)-1;
        subs_new(:,2,8) = subs(:,2)+1;
        
        %Consolidate into 2 dimensions
        subs_new_2 = cat(1,subs_new(:,:,1),subs_new(:,:,2));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,3));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,4));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,5));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,6));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,7));
        subs_new_2 = cat(1,subs_new_2,subs_new(:,:,8));
        
        subs_new = subs_new_2;
        
        %Throw away values outside range of image
        o1 = ones(size(subs_new(:,1))); %ones size row values
        o1(find(subs_new(:,1) > sz(1))) = nan; %too high
        o1(find(subs_new(:,1) == 0)) = nan; %too low
        
        o2 = ones(size(subs_new(:,2)));
        o2(find(subs_new(:,2) > sz(2))) = nan;
        o1(find(subs_new(:,2) == 0)) = nan;
        
        o = o1.*o2;
        delrows = find(isnan(o));
        subs_new(delrows,:) = [];
        
        %Keep only unique values
        inds_new = sub2ind(sz,subs_new(:,1),subs_new(:,2));
        [~,u_i,~] = unique(inds_new);
        inds_new = inds_new(u_i);
        
        %Throw out previously rejected points
        rej=ismember(inds_new,rejects);
        rejinds = find(rej);
        inds_new(rejinds) = [];
        
        %Throw out points already included
        inds_old = find(bimage);
        overlap = ismember(inds_new,inds_old);
        overlapinds = find(overlap);
        inds_new(overlapinds) = [];
        
        %Only keep points in cell
        inds_new = intersect(inds_new,cell_inds);
 
            if isempty(inds_new) == 1 %if there aren't points left, end
               grow = 0;
            else
                
                %Throw out points below cutoff value
                vals = image(inds_new);
                valdiff = vals-cutoff;
                out_low = find(valdiff<=0);
                
                %Throw out points ABOVE maximum cutoff value
                valdiff_high = vals-cutoff_high;
                out_high = find(valdiff_high>=0);
                
                out = [out_low ; out_high];
                
                
                if isempty(out) == 0
                    
                rejects = rejects(:);
                new_rejects = inds_new(out);
                new_rejects = new_rejects(:);
                rejects = cat(1,rejects,new_rejects);
                inds_new(out) = [];
                
                end
                
                %Set a new high cutoff if the number of iterations is 2 or
                %greater
                
                if iterations >= 2 
                    vals_new = image(inds_new);
                    cutoff_high = max(vals_new);
                    
                end
                
                if isempty(inds_new) == 1 %if there aren't points left, end
                    grow = 0;
                else

                    bimage(inds_new) = 1; %Set to 1 values in bimage
                    sz_inds = size(inds_new); %set up next loop :)
                    subs = zeros(sz_inds(1),2);
                    [subs(:,1),subs(:,2)] = ind2sub(sz,inds_new);
                    subs_new = zeros(sz_inds(1),2,8);
                    
                end
            end
              iterations = iterations+1;
     end
     
     clump_inds = find(bimage)';
     new_center = cluster_center;