function [bimage] = segment(image,point,parameter,exclude)
%Point in [X,Y] format, parameter -> 0 ,cutoff goes to global background, 
%parameter ->, cutoff goes to intensity of area directly around point
%Exclude- points that will be excluded from growth! Make nx1 vector, if you
%don't want to do it, just make it empty, 0x0 matrix

im = image;

%%%%%%%%%%%%% Image setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = double(im);
sz = size(im);
bimage = zeros(sz);

im = im - median(im(:));
im(im<0) = 0;


% parameter = 0.4;
im_gauss = imgaussfilt(im,4);
b_global = median(im_gauss(:));

% imagesc(im_gauss)
% [x,y] = ginput(1);
x = point(1);
y = point(2);
x = round(x);
y = round(y);
top=y-5;
bottom=y+5;
left=x-5;
right=x+5;
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

local_area=im_gauss(top:bottom,left:right);
i_cell = median(local_area(:));
cutoff = parameter*(i_cell-b_global)+b_global;
subs = [y,x];
bimage(subs(1),subs(2)) = 1;
grow = 1;
subs_new = zeros(1,2,8);
rejects = exclude;


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
        
        %Throw out points below cutoff value
        vals = im_gauss(inds_new);
        valdiff = vals-cutoff;
        out = find(valdiff<=0);
        
        if isempty(out) == 0
        new_rejects = inds_new(out);
        rejects = cat(1,rejects(:),new_rejects(:));
        inds_new(out) = [];
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
