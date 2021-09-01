function b1=compute_bas_v2(im,alp,pi,sig_2,nsamp,samps)

sz=size(im);
im=im(:);
sz1=numel(im);

ent=@(p) -1*(p.*log(p)+(1-p).*log(1-p));

t1=sigmoid(2*im/sig_2+log(pi)-log(1-pi));
t1p=sigmoid(2*im/sig_2-log(pi)+log(1-pi));




if isinf(nsamp(1))
    s11=ent(alp.*t1+(1-alp).*t1p);
else
    alp1=samps(1);
    s11=ent(alp1.*t1+(1-alp1).*t1p);
end
if isinf(nsamp(2))
    s21=alp.*ent(t1)+(1-alp).*ent(t1p);
else
    alp1=samps(2);
    s21=alp1.*ent(t1)+(1-alp1).*ent(t1p);
end




b1=reshape(s11-s21,sz);

end



