function filtered_image = med_filt_select(image,window_size,inds)

filtered_image = image;
sz = size(image);
nanim = nan(sz);
nanim(inds) = image(inds);

if rem(window_size,2) == 0
   hs = window_size/2; %half size
else
   hs = (window_size-1)/2;
end


for i  = 1:length(inds)
    
    [r,c] = ind2sub(sz,inds(i));
    top=r-hs;
    bottom=r+hs;
    left=c-hs;
    right=c+hs;
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
        
%     [rows,columns] = meshgrid(top:bottom,left:right);
%     area_inds = sub2ind(sz,rows,columns);
%     int_inds = intersect(inds,area_inds);
%     vals = image(int_inds);
%     filtered_image(inds(i)) = median(vals);

      local_area = nanim(top:bottom,left:right);
      vals = local_area(:);
      filtered_image(inds(i)) = median(vals,'omitnan');
    
end
