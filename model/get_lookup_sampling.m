sigs=logspace(-1,1,51);

parfor i=1:512
    i
    for j=1:16
        for k=1:51
            
            dmat=sign((dec2bin(i-1,9)-'0')-0.5);
            resp=j-1;
            lls(i,j,k)=ll_fixedsampling(sigs(k),resp,dmat,10000);
        end
    end
end