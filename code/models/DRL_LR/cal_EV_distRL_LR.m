function [opt1_expv,opt2_expv]=cal_EV_distRL_LR(alpha_pos_means,alpha_neg_means,vars1,vars2,opt1_events,opt2_events,start,nneuro)
%% inputs are base learning rate, betaa and betab for a beta distribution,events for opt1 and opt2,number of neurons simulated for the learnner
if nargin<8    nneuro=100; end
if nargin<7    start=0.5; end
%simulate a learner
betab_plr=(alpha_pos_means - vars1 + alpha_pos_means.*vars1 - alpha_pos_means.^2*2 + alpha_pos_means.^3)./vars1;
betaa_plr=-(alpha_pos_means.*(alpha_pos_means.^2 - alpha_pos_means + vars1))./vars1;

betab_nlr=(alpha_neg_means - vars2 + alpha_neg_means.*vars2 - alpha_neg_means.^2*2 + alpha_neg_means.^3)./vars2;
betaa_nlr=-(alpha_neg_means.*(alpha_neg_means.^2 - alpha_neg_means + vars2))./vars2;

alpha_pos=betarnd(betaa_plr,betab_plr,nneuro,1);
alpha_neg=betarnd(betaa_nlr,betab_nlr,nneuro,1);
%% learning two options
%learn first option
opt1_expv=ones(length(alpha_pos),size(opt1_events,1),size(opt1_events,2))*start;
for k=1:size(opt1_events,1)
    for i=1:size(opt1_events,2)-1
        for j=1:length(alpha_pos)
            if opt1_events(k,i)>opt1_expv(j,k,i)
             opt1_expv(j,k,i+1)=opt1_expv(j,k,i)+alpha_pos(j)*(opt1_events(k,i)-opt1_expv(j,k,i));
            else 
                opt1_expv(j,k,i+1)=opt1_expv(j,k,i)+alpha_neg(j)*(opt1_events(k,i)-opt1_expv(j,k,i));
            end
        end
    end
end
%learn the second option
opt2_expv=ones(length(alpha_pos),size(opt2_events,1))*start;
for k=1:size(opt2_events,1)
    for i=1:size(opt2_events,2)-1
        for j=1:length(alpha_pos)
            if opt2_events(k,i)>opt2_expv(j,k,i)
                opt2_expv(j,k,i+1)=opt2_expv(j,k,i)+alpha_pos(j)*(opt2_events(k,i)-opt2_expv(j,k,i));
            else
                opt2_expv(j,k,i+1)=opt2_expv(j,k,i)+alpha_neg(j)*(opt2_events(k,i)-opt2_expv(j,k,i));
            end
        end
    end
end