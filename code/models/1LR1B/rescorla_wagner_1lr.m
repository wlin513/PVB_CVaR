function [R]=rescorla_wagner_1lr(opt1_events,opt2_events,alpha,start)
%[r]=rescorla_wagner(Y,alpha,start)
% Y is column of wins and losses
% alpha is [reward_lr loss_lr], start is [start_reward start_loss]
%[r]=rescorla_wagner(Y,alpha)  (start is assumed to be 0.5
% Output is probability estimate

v1=zeros(size(opt1_events));
v1(:,1)=start;
v2=zeros(size(opt2_events));
v2(:,1)=start;
for b=1:size(opt1_events,1)
    for i=2:size(opt1_events,2)

        v1(b,i)=v1(b,i-1)+alpha*(opt1_events(b,i-1)-v1(b,i-1));
        
        v2(b,i)=v2(b,i-1)+alpha*(opt2_events(b,i-1)-v2(b,i-1));
    end
    
end

r=v1-v2;
R=reshape(r',size(r,1)*size(r,2),1);
