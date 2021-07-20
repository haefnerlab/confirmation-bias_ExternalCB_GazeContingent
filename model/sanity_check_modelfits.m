load('likelihood_4d_full');

sig_true=[10^(0),10^(0.5),10^(1)];
nsamps_true=[1,5,20];

nsamps=[1:2:100];
sigs=logspace(-1,2,51);

nreps=20;


for i=1:numel(sig_true)
    for j=1:numel(nsamps_true)
        for k=1:nreps
            num_trials=3e2;
            category_trial=sign(rand(num_trials,1)-0.5);
            dmat=category_trial.*sign(binornd(1,0.7,num_trials,9)-0.5);
            resp=simulate_model_ibs_v3([nsamps_true(j),sig_true(i)],dmat);
            ll_sig_ncog=get_like_sig_ncog(dmat,resp,lls1);
            ll1=ll_sig_ncog-logsumexp(ll_sig_ncog(:),1);
            ll1=exp(ll1);
            subplot(3,3,3*(i-1)+j);
            %         ll1=ll1';
            
            plot(nsamps,sum(ll1,2),'color',[1,0.75,0.75],'linewidth',2);
            
            ll1s{i,j}(k,:)=sum(ll1,2);
            
            hold on
            %         contour(nsamps,log10(sigs),ll1,3,'linewidth',2);
            
            %         xlim([1,100])
            %         ylim([-1,2])
            vline(nsamps_true(j),'m--');
            %         hline(log10(sig_true(i)),'m--');
            
            
            
            drawnow
        end
        
        
    end
end


for i=1:3
    for j=1:3
        subplot(3,3,3*(i-1)+j);
        plot(nsamps,mean(ll1s{i,j}),'r','linewidth',4);
        
        xlabel('num. of samples');
        ylabel('posterior probability');
        set(gca,'TickDir','out');
        set(gca,'linewidth',4);
        set(gcf, 'color', 'white')
        set(gca,'fontsize',20,'fontweight','bold')
        box off
        
        
    end
end


%%

nreps=5;
figure;
for i=1:numel(sig_true)
    for j=1:numel(nsamps_true)
        for k=1:nreps
        num_trials=1e4;
        category_trial=sign(rand(num_trials,1)-0.5);
        dmat=category_trial.*sign(binornd(1,0.7,num_trials,9)-0.5);
        resp=simulate_model_ibs_v3([nsamps_true(j),sig_true(i)],dmat);
        ll_sig_ncog=get_like_sig_ncog(dmat,resp,lls1);
        ll1=ll_sig_ncog-logsumexp(ll_sig_ncog(:),1);
        ll1=exp(ll1);
        subplot(3,3,3*(i-1)+j);
        %         ll1=ll1';
        
        plot(nsamps,sum(ll1,2),'color',[0.75,0.75,1],'linewidth',2);
        
        ll1s{i,j}(k,:)=sum(ll1,2);
        
        hold on
        %         contour(nsamps,log10(sigs),ll1,3,'linewidth',2);
        
        %         xlim([1,100])
        %         ylim([-1,2])
        vline(nsamps_true(j),'m--');
        %         hline(log10(sig_true(i)),'m--');
        
        
        
        drawnow
        
        xlabel('num. of samples');
        ylabel('posterior probability');
        set(gca,'TickDir','out');
        set(gca,'linewidth',4);
        set(gcf, 'color', 'white')
        set(gca,'fontsize',20,'fontweight','bold')
        box off
        end
    end
end


for i=1:3
    for j=1:3
        subplot(3,3,3*(i-1)+j);
        plot(nsamps,mean(ll1s{i,j}),'b','linewidth',4);
        
        xlabel('num. of samples');
        ylabel('posterior probability');
        set(gca,'TickDir','out');
        set(gca,'linewidth',4);
        set(gcf, 'color', 'white')
        set(gca,'fontsize',20,'fontweight','bold')
        box off
        
        
    end
end