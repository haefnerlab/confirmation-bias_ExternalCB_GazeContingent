function ll=ll_fixedsampling_v2(theta,dmat,nsamp)
tmp=99*ones(2.^size(dmat,2),1);
tmp1=ones(2.^size(dmat,2),nsamp);
for i=1:size(dmat,1)
    id=1+bin2dec(num2str(0.5*(1+dmat(i,:))));
    if tmp(id)==99
        resps(i,:)=simulate_model_ibs(theta,repmat(dmat(i,:),nsamp,1))';
        tmp(id)=1;
        tmp1(id,:)=resps(i,:);
    else
        resps(i,:)=tmp1(id,:);
    end
end
% parfor i=1:nsamp
%     resps(:,i)=simulate_model_ibs(theta,dmat);
% end

ll=zeros(1,16);
for j=1:size(resps,1)
    [q1,q2,q3]=histcounts(resps(j,:),[0:16]-0.5);
    pb=q1./nsamp;
    for k=1:16
    ll(k)=ll(k)+log(1e-4*(1/16)+(1-1e-4)*pb(k));
    end
end

end