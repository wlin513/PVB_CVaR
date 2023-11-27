%% load events
clear all
clc
%%
getfolders
addpath(genpath('utility\'))
addpath('models\PNPE')
addpath('models\1LR1b')
addpath('models\SI_model')
addpath('models\DRL')
addpath('models\CVaR')
addpath('models\DRL_LR')
addpath('models\PEIRS\')
addpath('models\Bayesian_CVaR\')
addpath('models\Bayesian_CVaR_eta_beta\')
addpath('models\Bayesian_CVaR_eta_priorvar\')
addpath('models\concave_UTIL\')
addpath('models\convex_UTIL\')
addpath('models\inverse_s_shape_UTIL\')
%%
load('events_for_nsim500.mat')
% load('events_v5.mat')
% 
% events.both_high=[BH_events(4,:);NH_events(3,:)];
% events.both_low=[BL_events(4,:);NL_events(3,:)];
% events.both_broad=[BH_events(1,:);BL_events(1,:)];
% events.both_narrow=[NH_events(1,:);NL_events(1,:)];
% 
% alltypes=fieldnames(events);
%%
%models={'Bayesian_CVaR_eta_beta','rw','PNPE','random','PEIRS','concave_UTIL','convex_UTIL','inverse_s_shape_UTIL'};%,'Bayesian_CVaR_eta_priorvar','PEIRS'
%models_names={'Bayesian-CVaR','1lr-RW','pos-neg-RW','random','PEIRS','concave-UTIL','convex-UTIL','inverse-s-shape-UTIL'};
model='Bayesian_CVaR';
mkdir([figdir,'sim_',model]);
tmpfigdir=[figdir,'sim_',model];
%
etabins=31;%bins for CVaR eta(positive bias/risk sensitivity)
betabins=30;
eta_space=-0.95:1.9/(etabins-1):0.95;
% log_beta_space=log(0.01):(log(20)-log(0.01))/(betabins-1):log(20);
% beta_space=exp(log_beta_space);
%beta_space=0.05:19.95/(betabins-1):20;
beta_space=5:45/(betabins-1):100;
%exp(log(0.05):(log(20)-log(0.05))/(betabins-1):log(20));

%
N=1000;%sample resolution
priora=1;%A value for prior beta distribution
priorb=1;%B value for prior beta distribution
updatevar=0.009;
step=1.01;
alpha=1;%Fixed Belif Model

nsim=500;
ntrial=30;
ifgenevents=0;
%%
    if ~ifgenevents
         load('events_for_nsim500.mat')
    end
%%
rt=([1:13]-1)./12.245+0.01;
f=figure;
for i=1:length(rt)
eventb=(rt(i) - updatevar + rt(i) .*updatevar - rt(i) .^2*2 + rt(i) .^3)./updatevar;
eventa=-(rt(i) .*(rt(i) .^2 - rt(i)  + updatevar))./updatevar;
fplot(@(x) betapdf(x,eventa,eventb),[0 1])
hold on;
ylim([0 5])
end
saveas(f,[figdir,'update_kernel.png'])
%%
for ns=1:nsim
    ns
    if ifgenevents
    %broad_high_deck
    v=0.04;%v=0.008;
    m=7.5/13;
    betab_bh=(m - v + m*v - 2*m^2 + m^3)/v;
    betaa_bh=-(m*(m^2 - m + v))/v;
    %broad_high_deck=12*betarnd(betaa_bh,betab_bh,[ntrial,1])+1;
    broad_high_deck=13*betarnd(betaa_bh,betab_bh,[ntrial,1])+0.5;
    broad_high_events(ns,:)=round(broad_high_deck);

    %narrow_high_deck
    v=0.01;
    m=7.5/13;
    betab_nh=(m - v + m*v - 2*m^2 + m^3)/v;
    betaa_nh=-(m*(m^2 - m + v))/v;
    narrow_high_deck=13*betarnd(betaa_nh,betab_nh,[ntrial,1])+0.5;
    narrow_high_events(ns,:)=round(narrow_high_deck);
    %broad_low_deck
    v=0.04;
    m=5.5/13;
    betab_bl=(m - v + m*v - 2*m^2 + m^3)/v;
    betaa_bl=-(m*(m^2 - m + v))/v;
    broad_low_deck=13*betarnd(betaa_bl,betab_bl,[ntrial,1])+0.5;
    broad_low_events(ns,:)=round(broad_low_deck);

    %narrow_low_deck
    v=0.01;
    m=5.5/13;
    betab_nl=(m - v + m*v - 2*m^2 + m^3)/v;
    betaa_nl=-(m*(m^2 - m + v))/v;
    narrow_low_deck=13*betarnd(betaa_nl,betab_nl,[ntrial,1])+0.5;
    narrow_low_events(ns,:)=round(narrow_low_deck);
    end
    events.both_high=[broad_high_events(ns,:);narrow_high_events(ns,:)];
    events.both_low=[broad_low_events(ns,:);narrow_low_events(ns,:)];
    events.both_broad=[broad_high_events(ns,:);broad_low_events(ns,:)];
    events.both_narrow=[narrow_high_events(ns,:);narrow_low_events(ns,:)];

    alltypes=fieldnames(events);

for i=1:length(alltypes)
    blktype=alltypes{i};
    data=(events.(blktype)'-1)/12.245+0.01;
    rrdist.(blktype)(:,:,1)=bayesian_update(data(:,1)',priora,priorb,alpha,updatevar,N);
    rrdist.(blktype)(:,:,2)=bayesian_update(data(:,2)',priora,priorb,alpha,updatevar,N);

   for j=1:length(eta_space)
       eta=eta_space(j);
       r.(blktype)(:,j,1)=readout_CVaR(rrdist.(blktype)(:,:,1),eta,step);
       r.(blktype)(:,j,2)=readout_CVaR(rrdist.(blktype)(:,:,2),eta,step);
   end
   for b=1:length(beta_space)
       beta=beta_space(b);
       pchoice.(blktype)(ns,:,:,b)=exp(beta*r.(blktype)(:,:,1))./(exp(beta*r.(blktype)(:,:,1))+exp(beta*r.(blktype)(:,:,2)));
   end
%end
%     %% figures
%     %events histogram
%     f1=figure;
%     histogram(events.(blktype)(1,:));hold on;histogram(events.(blktype)(2,:))
%     xlim([1 13])
%     xticks(1:13)
%     xlabel('values')
%     title(regexprep(blktype,'_',' '))
%     saveas(f1,[tmpfigdir,'\',blktype,'_event_hist.png'])
%     %all distribution estimates
%     f2=figure;
%     bcol=brewermap(length(data),"Blues");
%     rcol=brewermap(length(data),"Reds");
%     for  j=1:length(data)
%     plot(1:12/(N-1):13,squeeze(rrdist.(blktype)(j,:,1)),'Color',bcol(j,:))
%     hold on
%     plot(1:12/(N-1):13,squeeze(rrdist.(blktype)(j,:,2)),'Color',rcol(j,:))
%     end
%     xlabel('values')
%     xlim([1 13])
%     xticks(1:13)
%     ylabel('Probabilities')
%     title(regexprep(blktype,'_',' '))
%     saveas(f2,[tmpfigdir,'\',blktype,'.png'])
%     %gif showing the evolution of the distributions update
%     h = figure;
%     axis tight manual % this ensures that getframe() returns a consistent size
%     filename = [tmpfigdir,'\',blktype,'.gif'];
%     for j=1:length(data)
%     plot(1:12/(N-1):13,squeeze(rrdist.(blktype)(j,:,1)),'b')
%     hold on;
%     plot(1:12/(N-1):13,squeeze(rrdist.(blktype)(j,:,2)),'r')
%     hold off;
%     title(regexprep(blktype,'_',' '))
%     xlabel('values')
%     xlim([1 13])
%     xticks(1:13)
%     ylabel('Probabilities')
%     drawnow
%           % Capture the plot as an image 
%           frame = getframe(h); 
%           im = frame2im(frame); 
%           [imind,cm] = rgb2ind(im,256); 
%           % Write to the GIF File 
%           if j == 1 
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%           else 
%               imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%           end 
%     end
%     % trial by trial value update
%     f3=figure;
%     ax1=subplot(1,2,1);
%     plot((squeeze(r.(blktype)(:,:,1))-0.01)*12.245+1)
%     ylim([1 13])
%     yticks(1:13)
%     ylabel('values')
%     xlabel('trial n')
%     set(ax1,'Colororder',brewermap(etabins,"Blues"));
%     set(ax1,'Colormap',brewermap([],"Blues"));
%     c=colorbar;
%     c.Label.String='CVaR positive bias (eta)';
%     c.Label.Position=[0 0.5 0];
%     c.Label.Color='w';
%     c.Label.FontWeight='bold';
%     ax2=subplot(1,2,2);
%     plot((squeeze(r.(blktype)(:,:,2))-0.01)*12.245+1)
%     ylim([1 13])
%     yticks(1:13)
%     ylabel('values')
%     xlabel('trial n')
%     set(ax2,'Colororder',brewermap(etabins,"Reds"));
%     set(ax2,'Colormap',brewermap([],"Reds"));
%     c=colorbar;
%     c.Label.String='CVaR positive bias (eta)';
%     c.Label.Position=[0 0.5 0];
%     c.Label.Color='w';
%     c.Label.FontWeight='bold';
%     saveas(f3,[tmpfigdir,'\',blktype,'_trial_EVs.png'])
end
end
%%
f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(pchoice.(blktype)(:,:,:,end),1),2)),'LineWidth',5)
    hold on;
end
ylim([0.3 0.8])
%legend(regexprep(alltypes,'_','-'),'Location','southeastoutside','AutoUpdate','off')
xlabel('CVaR positive bias (eta)')
ylabel('probabilities of choosing the broader/higher option')
xticklabels(round([-0.95:1.9/6:0.95],2))
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')

[r,p]=corr(mean(optionchocies_all(sub_include_all,[5 7]),2),mean(optionchocies_all(sub_include_all,[6 8]),2),'type','Spearman')
[r,p]=corr(optionchocies_all(sub_include_all,7),optionchocies_all(sub_include_all,5),'type','Spearman')
figure;scatter(mean(optionchocies_all(sub_include_all,[5 7]),2),mean(optionchocies_all(sub_include_all,[6 8]),2))

figure;scatter(mean(optionchocies(sub_include,[5 6]),2),mean(optionchocies(sub_include,[7 8]),2))
[r,p]=corr(mean(optionchocies(sub_include,[5 6]),2),mean(optionchocies(sub_include,[7 8]),2),'type','Spearman')
[r,p]=corr(inv_logit(mean(optionchocies(sub_include,[5 6]),2)),inv_logit(mean(optionchocies(sub_include,[7 8]),2)))
[r,p]=corr(logit(mean(optionchocies(sub_include,[7 8]),2)),ques_table.SDS(sub_include),'type','Spearman')
[r,p]=corr(logit(optionchocies(sub_include,8)),ques_table.SDS(sub_include))
%%
% summary pro-v bias for each block type
f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(pchoice.(blktype),1),2),4)),'LineWidth',5)
    hold on;
end
ylim([0.3 0.8])
%legend(regexprep(alltypes,'_','-'),'Location','southeastoutside','AutoUpdate','off')
xlabel('CVaR positive bias (eta)')
ylabel('probabilities of choosing the broader/higher option')
xticklabels(round([-0.95:1.9/6:0.95],2))
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
saveas(f4,[tmpfigdir,'\provariance_for_each_blktype.png'])

f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(pchoice.(blktype),1),2),3)),'LineWidth',5)
    hold on;
end
%legend(regexprep(alltypes,'_','-'),'Location','southeastoutside','AutoUpdate','off')
xlabel('softmax temperature (beta)')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
xticklabels(round(0.01:19.99/6:20,2))
saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_beta.png'])
%%
save(['sim_results_',model,'.mat'])
