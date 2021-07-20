function ll_sig=get_like_sig(dmat,resp,lls)

ll_sig=zeros(size(lls,3),1);

for i=1:size(dmat,1)
    id_1=1+bin2dec(num2str(0.5*(1+dmat(i,:))));
    id_2=1+resp(i);
    
    ll_sig=ll_sig+squeeze(lls(id_1,id_2,:));
    
    
end


end