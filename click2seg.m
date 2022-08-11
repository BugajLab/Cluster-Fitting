%Click to segment :)

function [bimage,seed] = click2seg(image,parameter,exclude)
[x,y] = ginput(1);
seed = [x,y];
bimage = segment(image,seed,parameter,exclude);

