function [opt1_expv,opt2_expv]=cal_EV_distRL(alpha_base,betaa,betab,opt1_events,opt2_events,start,nneuro)
%% inputs are base learning rate, betaa and betab for a beta distribution,events for opt1 and opt2,number of neurons simulated for the learnner
if nargin<7    nneuro=100; end
if nargin<6    start=0.5; end
%simulate a learner
taus=betarnd(betaa,betab,nneuro,1);%%
% calculating alphas based on taus
%alpha_pos=2.*taus.*alpha_base;
%alpha_neg=2.*(1-taus)*alpha_base;

alpha_pos=2*taus*alpha_base;
alpha_pos(alpha_pos>=1)=0.999;
alpha_neg=2*(1-taus)*alpha_base;
alpha_neg(alpha_neg>=1)=0.999;
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