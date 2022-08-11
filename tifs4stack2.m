function [stack , filenames] = tifs4stack2(n)

[filenames,path] = uigetfile('*.tif','Select One or More Files','MultiSelect', 'on');

    
    stack = zeros(2048,2048,n);
    names = cell(n,1);
    for i = 1:n
    stack(:,:,i) = imread(filenames,i);
    names{i} = filenames;
    end
    
    filenames = names;
end