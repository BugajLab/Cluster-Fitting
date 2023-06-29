function [new_points] = growpoints(old_points,sz,iterations)

subs = zeros(length(old_points),2);
[subs(:,1),subs(:,2)] = ind2sub(sz,old_points);
bimage(old_points) = 1;
subs_new = zeros(length(old_points),2,8);


for i = 1:iterations
    
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
        
        %Throw out points already included
        inds_old = find(bimage);
        overlap = ismember(inds_new,inds_old);
        overlapinds = find(overlap);
        inds_new(overlapinds) = [];
  

                     bimage(inds_new) = 1; %Set to 1 values in bimage
                     sz_inds = size(inds_new); %set up next loop :)
                     subs = zeros(sz_inds(1),2);
                     [subs(:,1),subs(:,2)] = ind2sub(sz,inds_new);
                     subs_new = zeros(sz_inds(1),2,8);
end

new_points = find(bimage);