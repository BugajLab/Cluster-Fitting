function [stack , filenames] = tifs4stack()

[filenames,path] = uigetfile('*.tif','Select One or More Files','MultiSelect', 'on');





if isa(filenames,'cell') == 1
    
    filenames = filenames';
    szf = size(filenames);
    first_image = imread(filenames{1});
    sz = size(first_image);
    stack = zeros(sz(1),sz(2),szf(1));

    for i = 1:szf(1)
        stack(:,:,i) = imread(filenames{i});
    end

elseif isa(filenames,'char') == 1
    
    stack = imread(filenames);
    
end


