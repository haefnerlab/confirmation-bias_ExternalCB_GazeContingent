subjects = ...
    {
    'GCV3-subject01'; 'GCV3-subject02';...
    'GCV3-subject03'; 'GCV3-subject07';...
    'GCV3-subject08'; 'GCV3-subject11';...
    'GCV3-subject13'; 'GCV3-subject15';...
    'GCV3-subject16'; 'GCV3-subject17';...
    'GCV3-subject18'; 'GCV3-subject19';...
    'GCV3-subject20'; 'GCV3-subject21';...
    'GCV3-subject22'; 'GCV3-subject23';...
    };

load('likelihood_4d_full');
nsamps=[1:2:100];
sigs=logspace(-1,2,51);


%%
for kk=1:16
    kk
    data_sub = readSubjectData(subjects{kk});
    dmat1=data_sub.frame_categories;
    dmat=dmat1(:,[1,2,5,3,6,4,7,8,9]);
    resp=bin2dec(num2str([ones(numel(data_sub.choice),3),data_sub.choice(:)]));
    %%
    
    
    ll_sig_ncog=get_like_sig_ncog(dmat,resp,lls);
    ll1=ll_sig_ncog-logsumexp(ll_sig_ncog(:),1);
    ll1=exp(ll1);
    
    
    
    subplot(4,4,kk);
    %     plot(ns,ll_cog,'o-','linewidth',2)
    
    plot(log10(sigs),sum(ll1),'color',[1,0,0],'linewidth',2);
    
    xlabel('log_{10} \sigma_{peri}');
        ylabel('posterior probability');
        set(gca,'TickDir','out');
        set(gca,'linewidth',4);
        set(gcf, 'color', 'white')
        set(gca,'fontsize',20,'fontweight','bold')
        box off
    
    
    drawnow;
    
    
end