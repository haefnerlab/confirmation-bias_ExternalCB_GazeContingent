ns=[1:100];
load('sig_2_bfit_subs.mat');
sig_subs=unique(sig_2s_subs);


lls=zeros(512,16,100,6);
for i=1:512
    i
    
        
        for k=7:12
            k
            parfor l=1:numel(ns)
                
                dmat=sign((dec2bin(i-1,9)-'0')-0.5);
                resp=j-1;
                ss=sig_subs(k);
                lls(i,:,l,k)=ll_fixedsampling_v3([ns(l),sqrt(ss)],dmat,10000);
            end
        end
    
end