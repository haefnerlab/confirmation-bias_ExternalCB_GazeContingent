clear all;close all;clc;
figure(); hold on;
lims = [0 0 1920 1080];
center = [960 540];
num_peri = 3;
frames = 4;

frames = frames +1;
edge = (lims(4)-lims(2))/(sqrt(2)*frames);
pts = hexagonalGrid(lims,center,edge);
location = getLocation(lims,center,frames,num_peri);

colors = ['k','r','g','c','m'];
scatter(pts(:,1),pts(:,2),'b'); hold on;
for i=1:frames
f=location{i};
hold on;
scatter(f(:,1),f(:,2),colors(i))
pause;
% f2=location{2};
% hold on;
% pause;
% scatter(f2(:,1),f2(:,2),'g')
% f3 = location{3};
% hold on;
% pause;
% scatter(f3(:,1),f3(:,2),'k')
% f4 = location{4};
% hold on;
% pause;
% scatter(f4(:,1),f4(:,2),'c')
end