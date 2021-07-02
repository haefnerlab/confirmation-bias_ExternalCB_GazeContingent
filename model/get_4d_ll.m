sigs=logspace(-1,2,101);
nsamps=[1:100];
parfor i=1:512
    i
    
        for k=1:101
            k
            for l=1:100
            
            dmat=sign((dec2bin(i-1,9)-'0')-0.5);
            lls(i,:,l,k)=ll_fixedsampling_v3([nsamps(l),sigs(k)],dmat,10000);
            
            end
        end
    
end