sigs=logspace(-1,2,51);
nsamps=[81:2:100];
lls=zeros(512,16,numel(nsamps),numel(sigs));
for i=1:512
    i
    
        for k=1:numel(sigs)
            k
            parfor l=1:numel(nsamps)
            
            dmat=sign((dec2bin(i-1,9)-'0')-0.5);
            lls(i,:,l,k)=ll_fixedsampling_v3([nsamps(l),sigs(k)],dmat,10000);
            
            end
        end
    
end
    save(['/gpfs/fs1/home/sshivkum/ext_cb/confirmation-bias_ExternalCB_GazeContingent/model/likelihoods5'],'sigs','lls','nsamps');