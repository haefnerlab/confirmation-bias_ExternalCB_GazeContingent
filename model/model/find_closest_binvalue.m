function roundbin=find_closest_binvalue(bins,arr)
sz=size(arr);
bins=bins(:);
arr=arr(:)';

[~,closestIndex] = min(abs(bsxfun(@minus,bins, arr)));
roundbin=reshape(bins(closestIndex),sz);
end