function ll_cog=get_like_ncog(dmat,resp,lls)

ll_cog=zeros(size(lls,3),1);

for i=1:size(dmat,1)
    id_1=1+bin2dec(num2str(0.5*(1+dmat(i,:))));
    id_2=1+resp(i);
    
    ll_cog=ll_cog+squeeze(lls(id_1,id_2,:));
    
    
end


end