function ll_sig_ncog=get_like_sig_ncog(dmat,resp,lls)

ll_sig_ncog=zeros(size(lls,3),size(lls,4));

for i=1:size(dmat,1)
    id_1=1+bin2dec(num2str(0.5*(1+dmat(i,:))));
    id_2=1+resp(i);
    
    ll_sig_ncog=ll_sig_ncog+squeeze(lls(id_1,id_2,:,:));
    
    
end


end