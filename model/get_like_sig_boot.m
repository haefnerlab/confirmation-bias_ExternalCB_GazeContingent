function ll_sig_boot=get_like_sig_boot(dmat,resp,lls,nboot)


ll_sig_boot=zeros(nboot,size(lls,3));
dmat0=dmat;
resp0=resp;


parfor j=1:nboot
    j
    id=randi([1,size(dmat0,1)],[1,size(dmat0,1)]);
    
    resp=resp0(id);
    dmat=dmat0(id,:);
    
    
    ll_sig=zeros(size(lls,3),1);
    
    for i=1:size(dmat,1)
        id_1=1+bin2dec(num2str(0.5*(1+dmat(i,:))));
        id_2=1+resp(i);
        
        ll_sig=ll_sig+squeeze(lls(id_1,id_2,:));
        
        
    end
    
    ll_sig_boot(j,:)=ll_sig;
    
end


end