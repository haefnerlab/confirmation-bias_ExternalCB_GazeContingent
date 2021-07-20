function y = lognormpdf(X, mu, sigma)
sz=size(X);

X=X(:)';
mu=mu(:)';
sigs=(sigma(:))';
X=(X-mu)./sigs;
sigma=1;
mu=zeros(size(X));
d = size(X,1);
X = X-mu;
[U,p]= chol(sigma);
if p ~= 0
    error('ERROR: sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;

y=y-log(sigs);

y=reshape(y,sz);

end
