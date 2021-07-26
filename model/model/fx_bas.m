function y=fx_bas(eps1,alp,pi,sig_2,nsamp,sd,npts_leg)

sz=size(eps1);
eps1=eps1(:);
alp=alp(:);
alp=alp(:);
pi=pi(:);
sig_2=sig_2(:);
nsamp=nsamp(:);



[pts_leg, wts_leg] = GaussLegendre(npts_leg);


inds=cartprod(1:npts_leg,1:npts_leg);
ptt=[pts_leg(inds(:,1)),pts_leg(inds(:,2))];
wtt=[wts_leg(inds(:,1)),wts_leg(inds(:,2))];

im1=norminv(repmat(0.5+0.5*ptt(:,1),1,numel(eps1)),repmat(eps1(:)',npts_leg*npts_leg,1),sqrt(sig_2));
im2=norminv(repmat(0.5+0.5*ptt(:,2),1,numel(eps1)),repmat(-eps1(:)',npts_leg*npts_leg,1),sqrt(sig_2));
y1=(compute_bas_v3(im1,alp,pi,sig_2,nsamp,sd)>compute_bas_v3(im2,alp,pi,sig_2,nsamp,sd))+0.5*(compute_bas_v3(im1,alp,pi,sig_2,nsamp,sd)==compute_bas_v3(im2,alp,pi,sig_2,nsamp,sd));
y=0.5*0.5*sum(prod(wtt,2).*(y1));
y=reshape(y,sz);
end