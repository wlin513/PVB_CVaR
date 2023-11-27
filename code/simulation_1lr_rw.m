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
% load('events_v5.mat')
% 
% events.both_high=[BH_events(4,:);NH_events(3,:)];
% events.both_low=[BL_events(4,:);NL_events(3,:)];
% events.both_broad=[BH_events(1,:);BL_events(1,:)];
% events.both_narrow=[NH_events(1,:);NL_events(1,:)];
% 
% alltypes=fieldnames(events);
nsim=500;
ntrial=30;
%%
    %models={'Bayesian_CVaR_eta_beta','rw','PNPE','random','PEIRS','concave_UTIL','convex_UTIL','inverse_s_shape_UTIL'};%,'Bayesian_CVaR_eta_priorvar','PEIRS'
    %models_names={'Bayesian-CVaR','1lr-RW','pos-neg-RW','random','PEIRS','concave-UTIL','convex-UTIL','inverse-s-shape-UTIL'};
    %%
    model='1lr-RW';
    mkdir([figdir,'sim_',model]);
    tmpfigdir=[figdir,'sim_',model];
    %
    alphabins=30;%bins for CVaR eta(positive bias/risk sensitivity)
    betabins=30;
%     logit_alpha_space=logit(0.05):(logit(0.95)-logit(0.05))/(alphabins-1):logit(0.95);   
%     alpha_space=inv_logit(logit_alpha_space);
    alpha_space=0.05:0.9/(alphabins-1):0.95;
%     log_beta_space=log(0.01):(log(20)-log(0.01))/(betabins-1):log(20);
%     beta_space=exp(log_beta_space);
    beta_space=0.01:19.99/(betabins-1):20;
    start=0.5;
    
    ifgenevents=0;
    %%
    if ~ifgenevents
         load('events_for_nsim500.mat')
    end
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
    
    %%
    for i=1:length(alltypes)
        blktype=alltypes{i};
        data=(events.(blktype)'-1)/12.245+0.01;
        for j=1:length(alpha_space)
            alpha=alpha_space(j);
            r.(blktype)(:,j)=rescorla_wagner_1lr(data(:,1)',data(:,2)',alpha,start);
        end

       for b=1:length(beta_space)
           beta=beta_space(b);
           pchoice.(blktype)(ns,:,:,b)=1./(1+exp(-beta*r.(blktype)));
       end
    end
end
% summary pro-v bias for each block type
f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(pchoice.(blktype),1),2),4)),'LineWidth',5)
    hold on;
end
%legend(regexprep(alltypes,'_','-'),'Location','southeastoutside','AutoUpdate','off')
xticklabels(0.05:0.9/6:0.95)
xlabel('learning rate (alpha)')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
ylim([0.3 0.8])
saveas(f4,[tmpfigdir,'\provariance_for_each_blktype.png'])
%%
f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(pchoice.(blktype),1),2),3)),'LineWidth',5)
    hold on;
end
%legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('softmax temperature (beta)')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
xticklabels(round(0.01:19.99/6:20,2))
saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_beta.png'])
%%
% save('events_for_nsim500.mat','broad_high_events','broad_low_events','narrow_high_events','narrow_low_events')
save(['sim_results_',model,'.mat'])