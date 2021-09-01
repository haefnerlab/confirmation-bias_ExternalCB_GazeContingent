load('dump2');


%%
figure(1);


[x1,y1]=meshgrid(logspace(-1,log10(5),101),logspace(0,2,101));
vq=interp2(sigs,ns,acc,x1,y1,'linear');
vq1=interp2(sigs,ns,prop_confirmatory_saccade,x1,y1,'linear');
imagesc(vq)
set(gca,'Ydir','normal')
axis square
colormap gray
hold on
xticks([1:25:101])
yticks([1:25:101])
xticklabels(round(x1(1,[1:25:101]),1))
yticklabels(round(y1([1:25:101],1)/100,2))


se_acc=sqrt(subj_acc(:).*(1-subj_acc(:))./subj_ns(:));
se_prop=sqrt(subj_prop(:).*(1-subj_prop(:))./subj_ns(:));


[l,v,o]=isocontour(vq,subj_acc(2));
[~,id]=sort(v(:,1));
plot(v(id,2),v(id,1),'m-','linewidth',8);

if subj_acc(2)-se_acc(2)>min(vq(:))
    [l,v,o]=isocontour(vq,subj_acc(2)-se_acc(2));
    [~,id]=sort(v(:,1));
    plot(v(id,2),v(id,1),'m--','linewidth',4);
else
    plot(101*ones(size([1:101])),[1:101],'m--','linewidth',4);
end
if subj_acc(2)+se_acc(2)<max(vq(:))
    [l,v,o]=isocontour(vq,subj_acc(2)+se_acc(2));
    [~,id]=sort(v(:,1));
    plot(v(id,2),v(id,1),'m--','linewidth',4);
else
    plot(1*ones(size([1:101])),[1:101],'m--','linewidth',4);
end



[l,v,o]=isocontour(vq,subj_acc(1));
[~,id]=sort(v(:,1));
plot(v(id,2),v(id,1),'y-','linewidth',8);

if subj_acc(1)-se_acc(1)>min(vq(:))
    [l,v,o]=isocontour(vq,subj_acc(1)-se_acc(1));
    [~,id]=sort(v(:,1));
    plot(v(id,2),v(id,1),'y--','linewidth',4);
else
    plot(101*ones(size([1:101])),[1:101],'y--','linewidth',4);
end
if subj_acc(1)+se_acc(1)<max(vq(:))
    [l,v,o]=isocontour(vq,subj_acc(1)+se_acc(1));
    [~,id]=sort(v(:,1));
    plot(v(id,2),v(id,1),'y--','linewidth',4);
else
    
    plot(1*ones(size([1:101])),[1:101],'y--','linewidth',4);
end

pl(1)=plot(29,1,'r.','markersize',40);
pl(2)=plot(29,25,'g.','markersize',40);
pl(3)=plot(29,101,'b.','markersize',40);




c=colorbar;
c.Label.String ='accuracy';
xlabel({'peripheral uncertainty','( \sigma_{periphery} )'});
ylabel({'relative number of samples','( n_{cognitive}/n_{sensory} )'});
set(gca,'fontsize',40,'fontweight','normal')
set(gcf, 'color', 'white')

%%
figure(2);


[x1,y1]=meshgrid(logspace(-1,log10(5),101),logspace(0,2,101));
vq=interp2(sigs,ns,acc,x1,y1,'linear');
vq1=interp2(sigs,ns,prop_confirmatory_saccade,x1,y1,'linear');
imagesc(vq1)
set(gca,'Ydir','normal')
axis square
colormap gray
hold on
xticks([1:25:101])
yticks([1:25:101])
xticklabels(round(x1(1,[1:25:101]),1))
yticklabels(round(y1([1:25:101],1)/100,2))


se_acc=sqrt(subj_prop(:).*(1-subj_prop(:))./subj_ns(:));
se_prop=sqrt(subj_prop(:).*(1-subj_prop(:))./subj_ns(:));


[l,v,o]=isocontour(vq1,subj_prop(2));
[~,id]=sort(v(:,2));
plot(v(id,2),v(id,1),'m-','linewidth',8);
[l,v,o]=isocontour(vq1,max(min(vq1(:)),subj_prop(2)-se_prop(2)));
[~,id]=sort(v(:,2));
plot(v(id,2),v(id,1),'m--','linewidth',4);
[l,v,o]=isocontour(vq1,min(max(vq1(:)),subj_prop(2)+se_prop(2)));
[~,id]=sort(v(:,2));
plot(v(id,2),v(id,1),'m--','linewidth',4);


[l,v,o]=isocontour(vq1,subj_prop(1));
[~,id]=sort(v(:,2));
plot(v(id,2),v(id,1),'y-','linewidth',8);

if subj_prop(1)-se_acc(1)>min(vq1(:))
    [l,v,o]=isocontour(vq1,subj_prop(1)-se_prop(1));
    [~,id]=sort(v(:,2));
    plot(v(id,2),v(id,1),'y--','linewidth',4);
else
    plot(101*ones(size([1:101])),[1:101],'y--','linewidth',2);
end
if subj_prop(1)+se_acc(1)<max(vq1(:))
    [l,v,o]=isocontour(vq1,subj_prop(1)+se_prop(1));
    [~,id]=sort(v(:,2));
    plot(v(id,2),v(id,1),'y--','linewidth',4);
else
    
    plot(1*ones(size([1:101])),[1:101],'y--','linewidth',4);
end

pl(1)=plot(29,1,'r.','markersize',40);
pl(2)=plot(29,25,'g.','markersize',40);
pl(3)=plot(29,101,'b.','markersize',40);


c=colorbar;
c.Label.String ='prob. of confirmatory saccade';
xlabel({'peripheral uncertainty','( \sigma_{periphery})'});
ylabel({'relative number of samples','( n_{cognitive}/n_{sensory} )'});
set(gca,'fontsize',40,'fontweight','normal')
set(gcf, 'color', 'white')