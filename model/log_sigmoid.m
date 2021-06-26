function ly = log_sigmoid(x)
sz=size(x);
x=x(:);
ly = reshape(-logsumexp([zeros(size(x,1),1),-x],2),sz);
end