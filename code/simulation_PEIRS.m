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
%%
%models={'Bayesian_CVaR_eta_beta','rw','PNPE','random','PEIRS','concave_UTIL','convex_UTIL','inverse_s_shape_UTIL'};%,'Bayesian_CVaR_eta_priorvar','PEIRS'
%models_names={'Bayesian-CVaR','1lr-RW','pos-neg-RW','random','PEIRS','concave-UTIL','convex-UTIL','inverse-s-shape-UTIL'};
%%
model='PEIRS';
mkdir([figdir,'sim_',model]);
tmpfigdir=[figdir,'sim_',model];
%
alphaqbins=5;%bins for CVaR eta(positive bias/risk sensitivity)
alphasbins=5;
betabins=5;
omegabins=30;
spreadbins=5;
%logit_alpha_space=logit(0.05):(logit(0.95)-logit(0.05))/(alphaqbins-1):logit(0.95);   
%alpha_space_q=inv_logit(logit_alpha_space);

alpha_space_q=0.05:0.9/(alphasbins-1):0.95;


alpha_space_s=0.05:0.9/(alphasbins-1):0.95;

% if betabins>1
% log_beta_space=log(0.01):(log(20)-log(0.01))/(betabins-1):log(20);
% else 
%     log_beta_space=log(10);
% end
% beta_space=exp(log_beta_space);

beta_space=0.01:19.99/(betabins-1):20;


omega=-10:(10-(-10))/(omegabins-1):10;


s_log=log(0.01):(log(1)-log(0.01))/(spreadbins-1):log(1);

start=0.5;nsim=500;
ntrial=30;
    ifgenevents=0;
    %%
    if ~ifgenevents
         load('events_for_nsim500.mat')
    end
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
    for j=1:length(alpha_space_q)
        for k=1:length(alpha_space_q)
           for l=1:length(s_log)
             [r.(blktype)(:,j,k,l,1),r.(blktype)(:,j,k,l,2),s.(blktype)(:,j,k,l,1),s.(blktype)(:,j,k,l,2)]=PEIRS(data(:,1)',data(:,2)',alpha_space_q(j),alpha_space_q(k),start,exp(s_log(l)));
           end
        end
    end
pes=mean(r.(blktype),5)-0.5;%stimulus prediction errors
% replicate the value matrice to account for betas
for q=1:length(omega)
    q1(:,:,:,:,q)=r.(blktype)(:,:,:,:,1)+tanh(omega(q).*pes).*s.(blktype)(:,:,:,:,1);
    q2(:,:,:,:,q)=r.(blktype)(:,:,:,:,2)+tanh(omega(q).*pes).*s.(blktype)(:,:,:,:,2);
end
clear r s pes
val_l=q1-q2;
mmdl=repmat(val_l,[1,1,1,1,1,length(beta_space)]);

clear val_l q1 q2 

beta=repmat(permute(beta_space,[1 3 4 5 6 2]),[length(data) length(alpha_space_q) length(alpha_space_q) length(s_log) length(omega) 1]);

pchoice.(blktype)(ns,:,:,:,:,:,:)=1./(1+exp(-(beta.*mmdl)));

clear beta mmdl
end
end
%% summary pro-v bias for each block type
f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(mean(mean(mean(pchoice.(blktype),1),2),3),4),5),6)),'LineWidth',5)
    hold on;
end
legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('softmax temperature (beta)')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
%saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_beta.png'])

f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(mean(mean(mean(pchoice.(blktype),1),6),2),4),5),7)),'LineWidth',5)
    hold on;
end
legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('alpha Q')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
%saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_alpha_q.png'])

f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(mean(mean(mean(pchoice.(blktype),1),6),2),3),5),7)),'LineWidth',5)
    hold on;
end
legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('alpha S')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
%saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_alpha_s.png'])

f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(mean(mean(mean(pchoice.(blktype),1),6),2),3),4),7)),'LineWidth',5)
    hold on;
end
legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('start S')
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
%saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_starts.png'])

f4=figure;
for i=1:length(alltypes)
    blktype=alltypes{i};
    plot(squeeze(mean(mean(mean(mean(mean(mean(pchoice.(blktype),1),5),2),3),4),7)),'LineWidth',5)
    hold on;
end
%legend(regexprep(alltypes,'_','-'),'Location','northwestoutside','AutoUpdate','off')
xlabel('omega w')
xticklabels(round(-10:20/6:10,2))
ylabel('probabilities of choosing the broader/higher option')
title(strrep(model,'_','-'))
fplot(@(x) 0.5,[0 30],'--','Color','k')
ylim([0.3 0.8])
saveas(f4,[tmpfigdir,'\provariance_for_each_blktype_omega.png'])
%%
save(['sim_results_',model,'.mat'],'-v7.3')