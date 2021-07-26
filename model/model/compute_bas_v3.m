function b1=compute_bas_v3(im,alp,pi,sig_2,nsamp,sd)
rng(sd)
sz=size(im);
im=im(:);
alp=alp(:);
pi=pi(:);
sig_2=sig_2(:);
nsamp=nsamp(:);

sz1=numel(im);

% ent=@(p) -1*(p.*log(p)+(1-p).*log(1-p));
ent=@(log_p,log_1p) -1*(exp(log_p).*log_p+exp(log_1p).*log_1p);

log_t1=log_sigmoid(2*im./sig_2+log(pi)-log(1-pi));
log_t1p=log_sigmoid(2*im./sig_2-log(pi)+log(1-pi));

log_1t1=log_sigmoid(-2*im./sig_2-log(pi)+log(1-pi));
log_1t1p=log_sigmoid(-2*im./sig_2+log(pi)-log(1-pi));



if isinf(nsamp(1))
    s11=ent(logsumexp([log(alp)+log_t1,log(1-alp)+log_t1p],2),logsumexp([log(alp)+log_1t1,log(1-alp)+log_1t1p],2));
else
    rng(sd);
    alp1=binornd(repmat(nsamp(1),sz1,1),alp)/nsamp(1);
    s11=ent(logsumexp([log(alp1)+log_t1,log(1-alp1)+log_t1p],2),logsumexp([log(alp1)+log_1t1,log(1-alp1)+log_1t1p],2));
    
end
if isinf(nsamp(2))
    s21=alp.*ent(log_t1,log_1t1)+(1-alp).*ent(log_t1p,log_1t1p);
else
    alp1=binornd(repmat(nsamp(2),sz1,1),alp)/nsamp(2);
    s21=alp1.*ent(log_t1,log_1t1)+(1-alp1).*ent(log_t1p,log_1t1p);
end




b1=reshape(s11-s21,sz);

end



