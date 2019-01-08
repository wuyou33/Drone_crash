%% This script is used for print figures for publication.
function print_figures(handle,name,type,legend,resolution)

pos = get(gcf,'position');
size = get(gcf,'paperposition');

size(1:2) = zeros(1,2);
size(4) = size(3)/pos(3)*pos(4);
set(gcf,'paperposition',size);

print(handle,resolution,type,name);
end