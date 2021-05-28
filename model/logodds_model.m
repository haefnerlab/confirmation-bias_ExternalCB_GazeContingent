function ll=logodds_model(im,sig_2,pi)
ll=log((1-pi)+(pi).*exp(-2*im./sig_2))-log(pi+(1-pi).*exp(-2*im./sig_2));
end