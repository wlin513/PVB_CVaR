clear
clc
%%  set up event
load('../events.mat')
%1=Broad_high_event(BH_events);2=narrow_high_event(NH_events);
%3=Broad_low_event(BL_events);4=narrow_low_event(NL_events)
blocks_list=[1,2;3,4;1,4;1,3;2,3;2,4];
%blocks_list=[1,2;1,2;1,2;3,4;3,4;3,4;1,4;1,4;1,3;1,3;3,2;3,2;4,2;4,2];
[aBH,bBH]=find(blocks_list==1);
[aNH,bNH]=find(blocks_list==2);
[aBL,bBL]=find(blocks_list==3);
[aNL,bNL]=find(blocks_list==4);
for i=1:length(aBH)
    blocks_events(aBH(i),bBH(i),:)=BH_events(i,:);
    blocks_events(aNH(i),bNH(i),:)=NH_events(i,:);
    blocks_events(aBL(i),bBL(i),:)=BL_events(i,:);
    blocks_events(aNL(i),bNL(i),:)=NL_events(i,:);
end
%%
nsim=1000;
nabandon=5;
w=rand(nsim,1)*3.3-2.3;
w=exp(w);
%w=rand(nsim,1)*3;
% leak=rand(nsim,1)*0.5;
% noise=rand(nsim,1)*0.5;
% lapse=rand(nsim,1)*0.2;
leak=rand(nsim,1)*0.5;
% noise=rand(nsim,1)*3.3-2.3;
% noise=exp(noise);
noise=ones(nsim,1)*2;
lapse=rand(nsim,1)*0.2;
%
for i=0.1:0.1:3
    fplot(@(x) sigmf(x,[i,0]),[-3 3]);
    hold on
end
for i=1:nsim
    out=selective_gating(squeeze(blocks_events(:,1,:)),squeeze(blocks_events(:,2,:)),w(i),leak(i),noise(i),lapse(i),nabandon);
    % f1=figure;
    % plot(Y1);hold on;plot(Y2);hold on;plot(Y1-Y2)
    % legend('mid-broad','mid-narrow','diff')
    % title(['w=',num2str(w)])
    prob1(i,:,:)=out;

%mean(prob1)
%[h,p]=ttest(prob1-0.5)
end
%%
for i=1:size(prob1,1)
    for j=1:size(prob1,2)
        for k=1:size(prob1,3)
          choices(i,j,k)=randsample([1 0],1,1,[ prob1(i,j,k),1- prob1(i,j,k)]);       
        end
    end
end
positive_choices=choices(w>exp(mean(log(w))),:,:);
negative_choices=choices(w<exp(mean(log(w))),:,:);

blk_names={'block1:HB-HN','block2:LB-LN','block3:HB-LN','block4:HB-LB','block5:HN-LB','block6:HN-LN'};

f=figure('position',[63,122.5,1315,597.5]);
aa=[mean(mean(positive_choices,3),1);mean(mean(negative_choices,3),1)];
bb=[std(mean(positive_choices,3),0,1)/sqrt(size(positive_choices,1));std(mean(negative_choices,3),0,1)/sqrt(size(negative_choices,1))];
x=[0.85:1:5.85;1.15:1:6.15];
bar(aa'-0.5);
xticklabels(blk_names)
ylabel('choosing the narrow/making error ←      →correct choice/choosing the broad')
hold on
err=errorbar(x',aa'-0.5,bb');
err(1).LineStyle='none';
err(2).LineStyle='none';
err(1).LineWidth=2;
err(2).LineWidth=2;
err(1).Color='black';
err(2).Color='black';
legend('higher-inte-bias','lower-inte-bias','Location','SouthEast')
saveas(f,'barplot_choosing_broad.png')
%%
for i=1:6
    f=figure;

scatter(w(w>exp(mean(log(w)))),mean(choices(w>exp(mean(log(w))),i,:),3)-0.5,'b','filled')

hold on;
scatter(w(w<exp(mean(log(w)))),mean(choices(w<exp(mean(log(w))),i,:),3)-0.5,'r','filled')
% l=lsline;
% l(1).LineWidth=2;
% l(2).LineWidth=2;
% l(1).Color='r';
% l(2).Color='b';
xlabel('integration-bias')
title(blk_names(i))
saveas(f,['corr_integration_bias_block_',num2str(i),'.png'])
end
