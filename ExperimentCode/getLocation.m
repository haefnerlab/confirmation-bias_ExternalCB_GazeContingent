function location = getLocation(loc_limits,center,frames, num_peri)
edge_length = (loc_limits(4)-loc_limits(2))/(sqrt(2)*frames);
% pts = hexagonalGrid([480 30 1500 1050], [990 540], 170);
pts = hexagonalGrid([loc_limits(1) loc_limits(2) loc_limits(3) loc_limits(4)], [center(1) center(2)], edge_length);
location{1} = [center(1) center(2)];
for i=2:frames
    k = 1;
    for j=1:size(location{i-1},1)
        dist_norm = vecnorm((pts-location{i-1}(j,:))');
        b = dist_norm; b(b==min(b)) = [];
        nearest_pts = pts(round(dist_norm(:),2)==round(min(b),2),:);
        location{i}(k:k+num_peri-1,:) = nearest_pts(randperm(3,num_peri)',:);
        k = k+num_peri;
    end
    
end
end