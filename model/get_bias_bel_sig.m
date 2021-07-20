function bias=get_bias_bel_sig(bel,sig_2)
nsamp=1e5;
sig_2s=repmat(sig_2,nsamp,1);
im1=normrnd(1./sig_2s,sqrt(1./sig_2s));
im2=normrnd(-1./sig_2s,sqrt(1./sig_2s));
bas_diff=(compute_bas_v3(im1,bel,0.7,1,[Inf,Inf],randi(1e7))>compute_bas_v3(im2,bel,0.7,1,[Inf,Inf],randi(1e7)))+0.5*(compute_bas_v3(im1,bel,0.7,1,[Inf,Inf],randi(1e7))==compute_bas_v3(im2,bel,0.7,1,[Inf,Inf],randi(1e7)));
bias=mean(bas_diff);



end