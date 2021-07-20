function ll=ll_ibs(theta,resp,dmat,prop)
id=randperm(size(resp,1));
nn=max(1,floor(prop*size(resp,1)));
ll=0;
for i=1:nn
ll=ll+ibslike(@(a,b) simulate_model_ibs(a,b),theta,resp(id(i),:),dmat(id(i),:));
end
end