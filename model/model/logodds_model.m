function ll=logodds_model(im,sig_2,pi)
% ll=log((1-pi)+(pi).*exp(-2*im./sig_2))-log(pi+(1-pi).*exp(-2*im./sig_2));

if pi==0
    ll=(2*im./sig_2);
elseif pi==1
    ll=(-2*im./sig_2);
elseif im<-1e5
    ll=(log(pi)-log(1-pi))*ones(size(im));
elseif im>1e5
    ll=(log(1-pi)-log(pi))*ones(size(im));
else
    ll=log(1-pi)-log(pi)-log_sigmoid(2*im./sig_2-log(pi)+log(1-pi))+log_sigmoid(2*im./sig_2+log(pi)-log(1-pi));
end


end