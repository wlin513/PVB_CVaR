%% set the paths
clear all
close all
clc

parentdir='D:\OneDrive - University College London\2021_DistRL';
addpath(genpath([parentdir,'\code\models\']))
addpath(genpath([parentdir,'\code\utility\']))
figdir=[parentdir,'\final_figures\'];
%%
%dataset='combine_data'
%dataset='prolific';%'combine_data';%dataset='prolific';%%dataset='beijing';%dataset='fmri'
%dataset='discovery';
dataset='replication';
switch dataset
    case {'discovery'}
        datasetname='discovery sample';
        colors=[repmat([0.98,0.855,0.847],4,1);repmat([0.871,0.494,0.451],4,1)];
        load([parentdir,'\data\pro_variance_bias_prolific_discovery_sample.mat']);
    case {'replication'}
        datasetname='replication sample';
        colors=[repmat([0.68,0.85,0.92],4,1);repmat([0.415,0.686,0.902],4,1)];
        load([parentdir,'\data\pro_variance_bias_prolific_replication_sample.mat']);
    otherwise
        disp('That one does not exist!')
end

%%
start=0.5;
abandontn=0;
ntrials=30;
%%
ques=struct2table(data);
[optionchocies,sub_include]=calculate_prov(data,ntrials);
%%
tic
for i=1:size(data,2)
     disp(i) 
      choices=transpose(reshape(data(i).choice,ntrials,length(data(i).blockN)));
      %blk=find(data(i).blockN==5|data(i).blockN==6|data(i).blockN==7|data(i).blockN==8);
      blk=[find(data(i).blockN==5),find(data(i).blockN==6),find(data(i).blockN==7),find(data(i).blockN==8)];

      choice=choices(blk,1:ntrials);
      resp_made=true(size(choice));

  
    opt1_events=transpose(reshape(data(i).opt1_out,ntrials,length(data(i).blockN)));
    opt1_events=opt1_events(blk,1:ntrials);
    opt1_events=(opt1_events-1)/12.245+0.01;

    opt2_events=transpose(reshape(data(i).opt2_out,ntrials,length(data(i).blockN)));
    opt2_events=opt2_events(blk,1:ntrials);
    opt2_events=(opt2_events-1)/12.245+0.01;
%     estimates_Bayesian_CVaR_eta_beta(i)=fit_linked_Bayesian_CVaR_eta_beta(opt1_events,opt2_events,choice,abandontn,resp_made,1,30,30,0);
%      estimates_1LR1B(i)=fit_linked_1LR1B(opt1_events,opt2_events,choice,start,abandontn);
%      estimates_PNPE(i)=fit_linked_PNPE(opt1_events,opt2_events,choice,start,abandontn);
%      estimates_PEIRS(i)=fit_linked_PEIRS(opt1_events,opt2_events,choice,start,abandontn);
%      estimates_concave_UTIL(i)=fit_linked_concave_UTIL(opt1_events,opt2_events,choice,start,abandontn,resp_made,30,30,30,0);
     estimates_convex_UTIL(i)=fit_linked_convex_UTIL(opt1_events,opt2_events,choice,start,abandontn,resp_made,30,30,30,0);
%      estimates_inv_s_shape_UTIL(i)=fit_linked_inverse_s_shape_UTIL(opt1_events,opt2_events,choice,start,abandontn,resp_made,30,30,30,0);
%     resp_made=resp_made(:,abandontn+1:end);

    random(i).BIC=2*-sum(sum(log(0.5*resp_made+1e-16)));%change this if there is actually no response trials
end
toc
ee.Bayesian_CVaR_eta_beta=struct2table(estimates_Bayesian_CVaR_eta_beta);
%ee.Bayesian_CVaR_eta_priorvar=struct2table(estimates_Bayesian_CVaR_eta_priorvar);
ee.rw=struct2table(estimates_1LR1B);
ee.PNPE=struct2table(estimates_PNPE);
ee.concave_UTIL=struct2table(estimates_concave_UTIL);
ee.convex_UTIL=struct2table(estimates_convex_UTIL);
ee.inverse_s_shape_UTIL=struct2table(estimates_inv_s_shape_UTIL);
ee.random=struct2table(random);
ee.PEIRS=struct2table(estimates_PEIRS);

pos_bias=ee.PNPE.mean_alpha_pos./(ee.PNPE.mean_alpha_pos+ee.PNPE.mean_alpha_neg);
pos_bias_logit=logit(pos_bias);
clear PEIRS CVAR PNPE RW
for i=1:size(estimates_PEIRS,2)
    PEIRS(i,:)=mean(estimates_PEIRS(i).prob_ch_left,2);
    CVAR(i,:)=mean(estimates_Bayesian_CVaR_eta_beta(i).prob_ch_left,2);
    PNPE(i,:)=mean(reshape(estimates_PNPE(i).prob_ch_left,30,4),1);
    RW(i,:)=mean(reshape(estimates_1LR1B(i).prob_ch_left,30,4),1);
    convex(i,:)=mean(reshape(estimates_convex_UTIL(i).prob_ch_left,30,4),1);
end
%%
save(['model_results_',dataset,'.mat'],  '-v7.3')
%%
marg_all=[];
for i=1:size(estimates_Bayesian_CVaR_eta_beta,2)
    if sub_include(i)==1
    marg_all=[marg_all,estimates_1LR1B.marg_beta];
    end
end
figure;plot(marg_all);hold on;plot(mean(marg_all,2),'LineWidth',5)
%%
%figure;plot(smooth(squeeze(mean(mean(opt1chocies_all(sub_include_all,[5 7],:)))),6));hold on;plot(smooth(squeeze(mean(mean(opt1chocies_all(sub_include_all,[6 8],:)))),6))
sf1=figure;
set(gcf,'Position',[300 300 600 370])
subplot(1,2,1)
plot(smooth(squeeze(mean(opt1chocies_all(sub_include_all,5,:))),6),'LineWidth',5);hold on;plot(smooth(squeeze(mean(opt1chocies_all(sub_include_all,6,:))),6),'LineWidth',5)
hold on
fplot(@(x) 0.5+0*x,[0 30],'--','Color','k')
legend('BHNH','BLNL','Location','northwest')
legend boxoff
xlabel('trials')
ylabel('percentages of participants choosing the broader option')
subplot(1,2,2)
plot(smooth(squeeze(mean(opt1chocies_all(sub_include_all,7,:))),6),'LineWidth',5);hold on;plot(smooth(squeeze(mean(opt1chocies_all(sub_include_all,8,:))),6),'LineWidth',5)
hold on
fplot(@(x) 0.5+0*x,[0 30],'--','Color','k')
legend('BiHNH','BiLNL','Location','northwest')
legend boxoff
xlabel('trials')
ylabel('percentages of participants choosing the broader option')
exportgraphics(sf1,[figdir,'trialbytrial_PVB.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(sf1,[figdir,'trialbytrial_PVB.png'],'BackgroundColor','none','Resolution',300)
%% figure 1 
f1=figure;
set(gcf,'Position',[300 300 370 250])
for i=1:8
b=bar(i,mean(optionchocies(sub_include,i),1));

set(b,'FaceColor',colors(i,:),'EdgeColor',colors(i,:))
 hold on
     dotPlot_xtr(optionchocies(sub_include,i),i,colors(9-i,:),0.05)
end
ylabel('precentages of choosing option a','FontSize',10);
xticks([1:8])
xticklabels({'BHBL','NHNL','BHNL','NHBL','BHNH','BLNL','BiHNH','BiLNL'});
title(datasetname,'FontSize',10)
ax=gca;
ax.XAxis.FontSize=5;
ax.YAxis.FontSize=8;
hold on
err=std(optionchocies(sub_include,:))/sqrt(sum(sub_include));
er=errorbar(mean(optionchocies(sub_include,:),1),err,'Color','black');
er.LineStyle='none';
hold on
fplot(@(x) 0.5+0*x,[-0.2 9.2],'--','Color','k')
exportgraphics(f1,[figdir,dataset,'bar_perf_block.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f1,[figdir,dataset,'bar_perf_block.png'],'BackgroundColor','none','Resolution',300)

f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(ques.RRS(sub_include),mean(optionchocies(sub_include,5:8),2),1,1,'ruminative response scales','provariance-bias',colors(end,:));
title(datasetname,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_RRS_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_RRS_pro_v.png'],'BackgroundColor','none','Resolution',300)

%% for combined data
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(option_include_all&accuracy_include_all),mean(optionchocies_all(option_include_all&accuracy_include_all,5:8),2),1,1,'ruminative response scales','provariance-bias',colors(end,:));
title(datasetname,'FontSize',10)
dataset='combine'
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(sub_include_all),mean(optionchocies_all(sub_include_all,5),2),1,1,'ruminative response scales','PVB-BHNH',colors(end,:));
title(dataset,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_BHNH_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_BHNH_pro_v.png'],'BackgroundColor','none','Resolution',300)
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(sub_include_all),mean(optionchocies_all(sub_include_all,6),2),1,1,'ruminative response scales','PVB-BLNL',colors(end,:));
title(dataset,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_BLNL_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_BLNL_pro_v.png'],'BackgroundColor','none','Resolution',300)
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(sub_include_all),mean(optionchocies_all(sub_include_all,7),2),1,1,'ruminative response scales','PVB-BiHNH',colors(end,:));
title(dataset,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_BiHNH_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_BiHNH_pro_v.png'],'BackgroundColor','none','Resolution',300)
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(sub_include_all),mean(optionchocies_all(sub_include_all,8),2),1,1,'ruminative response scales','PVB-BiLNL',colors(end,:));
title(dataset,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_BiLNL_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_BiLNL_pro_v.png'],'BackgroundColor','none','Resolution',300)
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(RRS(option_include_all&accuracy_include_all),mean(optionchocies_all(option_include_all&accuracy_include_all,5:8),2),1,1,'ruminative response scales','provariance-bias',colors(end,:));
title(datasetname,'FontSize',10)

exportgraphics(f3,[figdir,dataset,'_RRS_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f3,[figdir,dataset,'_RRS_pro_v.png'],'BackgroundColor','none','Resolution',300)
%%
[r,p]=corr(mean(optionchocies(sub_include,5:8),2),logit(ee.rw.mean_alpha(sub_include)))
[r,p]=corr(mean(optionchocies(sub_include,5:8),2),ee.PEIRS.mean_omega(sub_include))

[r,p]=corr(ques.RRS(sub_include),logit(ee.rw.mean_alpha(sub_include)))
[r,p]=corr(ques.RRS(sub_include),ee.PEIRS.mean_omega(sub_include))
%% figure 5
f2=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(mean(optionchocies(sub_include,5:8),2),logit((ee.Bayesian_CVaR_eta_beta.mean_eta(sub_include)+1)/2),1,1,'provariance-bias','risk aversion\leftarrow CVaR\rightarrow risk seeking',colors(end,:));
title(datasetname,'FontSize',10)
exportgraphics(f2,[figdir,dataset,'_eta_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
exportgraphics(f2,[figdir,dataset,'_eta_pro_v.png'],'BackgroundColor','none','Resolution',300)
f2=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(mean(optionchocies(sub_include,5:8),2),pos_bias_logit(sub_include),1,1,'provariance-bias','pos-bias-lr',colors(end,:));
title(datasetname,'FontSize',10)
exportgraphics(f2,[figdir,dataset,'_pos_bias_lr_pro_v.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f2,[figdir,dataset,'_pos_bias_lr_pro_v.png'])
f3=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(ques.RRS(sub_include),pos_bias_logit(sub_include),1,1,'ruminative response scales','pos-bias-lr',colors(end,:));
title(datasetname,'FontSize',10)
exportgraphics(f3,[figdir,dataset,'_RRS_pos_lr.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f3,[figdir,dataset,'_RRS_por_lr.png'])
f4=figure;
set(gcf,'Position',[300 300 250 250])
plotcorr(ques.RRS(sub_include),logit((ee.Bayesian_CVaR_eta_beta.mean_eta(sub_include)+1)/2),1,1,'ruminative response scales','risk aversion\leftarrow CVaR\rightarrow risk seeking',colors(end,:));
title(datasetname,'FontSize',10)
exportgraphics(f4,[figdir,dataset,'_RRS_eta.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f4,[figdir,dataset,'_RRS_eta.png'])

f5=figure;
set(gcf,'Position',[300 300 370 250])
models={'Bayesian_CVaR_eta_beta','rw','PNPE','convex_UTIL','PEIRS','concave_UTIL','random','inverse_s_shape_UTIL'};%,'Bayesian_CVaR_eta_priorvar','PEIRS'
models_names={'Bayesian-CVaR','1lr-RW','pos-neg-RW','convex-UTIL','PEIRS','concave-UTIL','random','inverse-s-shape-UTIL'};
%models={'Bayesian_CVaR_eta_beta','rw','PNPE','random','PEIRS'};%,'Bayesian_CVaR_eta_priorvar','PEIRS'
%models_names={'Bayesian-CVaR','1lr-RW','pos-neg-RW','random','PEIRS'};
BICs=[];
for i=1:length(models)
BICs=[BICs;sum(ee.(models{i}).BIC(sub_include))];
end
barh(BICs-min(BICs),'FaceColor',colors(end,:),'EdgeColor',colors(end,:))
title(datasetname,'FontSize',10)
ax=gca;
ax.XAxis.FontSize=5;
ax.YAxis.FontSize=8;
yticklabels(models_names)
xlabel('∆BICs','FontSize',8)
exportgraphics(f5,[figdir,dataset,'_delta_BICs.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f5,[figdir,dataset,'_delta_BICs.png'])
%% plot pro-v from simulations using best-fited parameter for each model
f5=figure;
set(gcf,'Position',[300 300 370 370])
drw=mean(RW(sub_include,[1 3]),2)-mean(RW(sub_include,[2 4]),2);
dpnpe=mean(PNPE(sub_include,[1 3]),2)-mean(PNPE(sub_include,[2 4]),2);
dpeirs=mean(PEIRS(sub_include,[1 3]),2)-mean(PEIRS(sub_include,[2 4]),2);
dcvar=mean(CVAR(sub_include,[1 3]),2)-mean(CVAR(sub_include,[2 4]),2);
dconvex=mean(convex(sub_include,[1 3]),2)-mean(convex(sub_include,[2 4]),2);
ddata=mean(optionchocies(sub_include,[5 7]),2)-mean(optionchocies(sub_include,[6 8]),2);
bar(mean([dpnpe,dconvex,drw,dpeirs,dcvar,ddata]),'FaceColor',colors(end,:),'EdgeColor',colors(end,:))
hold on;
errorbar(mean([dpnpe,dconvex,drw,dpeirs,dcvar,ddata]),std([dpnpe,dconvex,drw,dpeirs,dcvar,ddata])/sqrt(sum(sub_include)),'LineStyle','none')
xticklabels({'pos-neg-RW','convex-UTIL','1lr-RW','PEIRS','Bayesian-CVaR','data'})
xtickangle(45)
ylim([-0.05 0.4])
box off
ylabel('P(risky|both-high) - P(risky|both-low)')
exportgraphics(f5,[figdir,dataset,'_sim_delta_equalblks.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f5,[figdir,dataset,'_sim_delta_equalblk.png'])
figure;bar(mean(CVAR))
figure;bar(mean(PEIRS))
figure;bar(mean(PNPE))
[h,p,~,stat]=ttest2(dcvar,ddata)
[h,p,~,stat]=ttest2(dpeirs,ddata)
[h,p,~,stat]=ttest2(drw,ddata)
[h,p,~,stat]=ttest2(dconvex,ddata)
[h,p,~,stat]=ttest2(dpnpe,ddata)

figure;scatter(mean(PNPE,2)-0.5,ee.PNPE.mean_alpha_pos./(ee.PNPE.mean_alpha_pos+ee.PNPE.mean_alpha_neg))
figure;scatter(mean(CVAR,2)-0.5,ee.Bayesian_CVaR_eta_beta.mean_eta)
figure;scatter(CVAR(:,4),ee.Bayesian_CVaR_eta_beta.mean_eta)
figure;scatter(log(ee.Bayesian_CVaR_eta_beta.mean_beta),ee.Bayesian_CVaR_eta_beta.mean_eta)

[h,p]=ttest(PEIRS(sub_include,1)-PEIRS(sub_include,2))
[h,p,stat]=ttest(PNPE(sub_include,1)-PNPE(sub_include,2))
[h,p,stat]=ttest(CVAR(sub_include,1)-CVAR(sub_include,2))
[h,p,stat]=ttest(RW(:,1)-RW(:,2))

figure;scatter(ee.Bayesian_CVaR_eta_beta.mean_eta(sub_include),CVAR(sub_include,1)-CVAR(sub_include,2))
figure;scatter(ee.Bayesian_CVaR_eta_beta.mean_eta(sub_include),CVAR(sub_include,1)-CVAR(sub_include,2))



%% save data as csv file for spss
tt=optionchocies(sub_include,:);
T=array2table(tt);
T.Properties.VariableNames={'BHBL','NHNL','BHNL','NHBL','BHNH','BLNL','BiHNH','BiLNL'};
T.PVB=mean(optionchocies(sub_include,5:8),2);
T.DAS=ques.DAS(sub_include);
T.RRS=ques.RRS(sub_include);
T.SDS=ques.SDS(sub_include);
T.tSTAI=ques.tSTAI(sub_include);
T.IUS=ques.IUS(sub_include);
writetable(T,[dataset,'_pro_variance_bias_each_blk.csv'])
%%
load([parentdir,'\data\pro_variance_bias_prolific_discovery_sample.mat']);
data1=data;
load([parentdir,'\data\pro_variance_bias_prolific_replication_sample.mat']);
data2=data;
ques1=struct2table(data1);
ques2=struct2table(data2);
[optionchocies1,sub_include1]=calculate_prov(data1,ntrials);
[optionchocies2,sub_include2]=calculate_prov(data2,ntrials);
colors=[0.871,0.494,0.451;0.415,0.686,0.902];
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques1.RRS(sub_include1),'FaceColor',colors(1,:));hold on;histogram(ques2.RRS(sub_include2),'FaceColor',colors(2,:));
title('Ruminaton Response Scale')
leg=legend({'discovery sample','replication sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'histogram_RRS.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'histogram_RRS.png'])
[h,p]=ttest2(ques1.RRS(sub_include1),ques2.RRS(sub_include2))
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques1.SDS(sub_include1),'FaceColor',colors(1,:),'BinWidth',10);hold on;histogram(ques2.SDS(sub_include2),'FaceColor',colors(2,:),'BinWidth',10);
[h,p]=ttest2(ques1.SDS(sub_include1),ques2.SDS(sub_include2))
title('Self-Rating Depression Scale')
leg=legend({'discovery sample','replication sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'histogram_SDS.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'histogram_SDS.png'])
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques1.tSTAI(sub_include1),'FaceColor',colors(1,:));hold on;histogram(ques2.tSTAI(sub_include2),'FaceColor',colors(2,:));
[h,p]=ttest2(ques1.tSTAI(sub_include1),ques2.tSTAI(sub_include2))
title('Trait Anxiety')
leg=legend({'discovery sample','replication sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'histogram_tSTAI.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'histogram_tSTAI.png'])
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques1.IUS(sub_include1),'FaceColor',colors(1,:));hold on;histogram(ques2.IUS(sub_include2),'FaceColor',colors(2,:));
title('Intolerance of Uncertainty Scale')
leg=legend({'discovery sample','replication sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'histogram_IUS.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'histogram_IUS.png'])
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques1.DAS(sub_include1),'FaceColor',colors(1,:));
title('Dysfunctional Attitude Scale form A')
leg=legend({'discovery sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'discovery','_histogram_DAS.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'discovery','histogram_DAS.png'])
f8=figure;
set(gcf,'Position',[300 300 350 250])
histogram(ques2.HAPPI(sub_include2),'FaceColor',colors(2,:));
title('Hypomanic Attitudes and Positive Predictions Inventory','FontSize',8)
leg=legend({'discovery sample','replication sample'},'Box','off','FontSize',6);
leg.ItemTokenSize=[5,5];
exportgraphics(f8,[figdir,'replication','_histogram_HAPPI.eps'],'BackgroundColor','none','ContentType','vector')
saveas(f8,[figdir,'replication','_histogram_HAPPI.png'])
[h,p]=ttest2(ques1.IUS(sub_include1),ques2.IUS(sub_include2))

qq=[ques1.RRS;ques2.RRS];
prov=[mean(optionchocies1(:,6),2);mean(optionchocies2(:,6),2)];
sinc=[sub_include1;sub_include2];
[r,p]=corr(qq(sinc),prov(sinc))