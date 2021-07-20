%%
num_reps=3;


sig_2=@(alp,sige_2) (1./(norminv(alp).^2))-sige_2;

alps=linspace(0.5,1,11);
alps([1,end])=[0.51,0.99];
alps=0.8;
num_sacs=[0:15];

for q4=1:numel(num_sacs)
    num_sac=num_sacs(q4);
    
    for q3=1:numel(alps)
        
        sig_peri=sqrt(sig_2(alps(q3),0.11.^2));
        
        
        for q1=1:num_reps
            clc
            q4
            q3
            q1
            num_trials=1e5;
            category_trial=sign(rand(num_trials,1)-0.5);
            num_images=3+2*num_sac;
            p_match=0.7;
            sig_fov=sqrt(0.01.^2+0.11.^2);
            %     sig_peri=sqrt((3).^2+0.11.^2);
            frame_categories=repmat(category_trial,1,num_images).*sign(binornd(1,p_match,num_trials,num_images)-0.5);
            im_fov=(frame_categories+normrnd(0,sig_fov,num_trials,num_images));
            im_peri=(frame_categories+normrnd(0,sig_peri,num_trials,num_images));
            
            gamma=(sig_peri.^2)./((sig_peri.^2)+(sig_fov.^2));
            lo_fov=-1*logodds_model(im_fov,sig_fov.^2,p_match);
            lo_fov1=-1*logodds_model(im_fov*gamma+im_peri*(1-gamma),(sig_fov.^2)*gamma,p_match);
            lo_peri=-1*logodds_model(im_peri,sig_peri.^2,p_match);
            
            %%
            for q2=[1,6,11]%[1:11]
                %         q2
                p_bias=0.1*(q2-1);
                
                lo=lo_fov(:,1);
                %         locs=[2*ones(num_trials,1),4*ones(num_trials,1),6*ones(num_trials,1),8*ones(num_trials,1),10*ones(num_trials,1),12*ones(num_trials,1),14*ones(num_trials,1),16*ones(num_trials,1),18*ones(num_trials,1),20*ones(num_trials,1)];
                locs=[];
                for i1=1:num_sac
                    locs=[locs,2*i1*ones(num_trials,1)];
                end
                for i=1:num_sac
                    %             lo=lo+sum(lo_peri(:,[2+2*(i-1),3+2*(i-1)]),2);
                    l1=2*(i-1)+2;
                    l2=2*(i-1)+3;
                    swap=rand(num_trials,1)>p_bias;
                    
                    ids2=(sign(im_peri(:,l1))==sign(im_peri(:,l2))) | (0==sign(lo));
                    ids2=ids2 & (rand(num_trials,1)<0.5);
                    locs(ids2,i)=l2;
                    
                    
                    %             ids=(~ids2)&(im_peri(:,l2).*sign(lo)>im_peri(:,l1).*sign(lo));
                    
                    ids=(~ids2)&(sign(im_peri(:,l2))==sign(lo));
                    locs(ids,i)=l2;
                    
                    
                    locs((~ids2)&swap & locs(:,i)==l1,i)=100;
                    locs((~ids2)&swap & locs(:,i)==l2,i)=l1;
                    locs((~ids2)&swap & locs(:,i)==100,i)=l2;
                    lo1=zeros(num_trials,1);
                    lo1(locs(:,i)==l1)=lo_fov1(locs(:,i)==l1,l1);
                    lo1(locs(:,i)==l2)=lo_fov1(locs(:,i)==l2,l2);
                    
                    lo2=zeros(num_trials,1);
                    lo2(locs(:,i)==l1)=lo_peri(locs(:,i)==l1,l2);%+lo_peri(locs(:,i)==l1,l1);
                    lo2(locs(:,i)==l2)=lo_peri(locs(:,i)==l2,l1);%+lo_peri(locs(:,i)==l2,l2);
                    
                    
                    lo=lo+lo1+lo2;
                    
                end
                
                lo=lo+sum(lo_peri(:,[end-1,end]),2);
                acc(q1,q2,q3,q4)=mean(sign(lo)==category_trial);
%                 idd=(sign(im_peri(:,2))~=sign(im_peri(:,3)));%  & (sign(im_fov(:,1))==sign(im_peri(:,2))) & (sign(im_fov(:,3))==sign(im_peri(:,3)));% & (sign(im_fov(:,3))==sign(im_peri(:,3)));% & (sign(im_fov(:,3))==sign(im_peri(:,3)));
%                 acc1(q1,q2,q3,q4)=mean(sign(lo(idd))==category_trial(idd));
%                 acc2(q1,q2,q3,q4)=mean(sign(lo(~idd))==category_trial(~idd));
            end
        end
        
        
    end
end