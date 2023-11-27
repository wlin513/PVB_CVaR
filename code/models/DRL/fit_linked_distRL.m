function [out]=fit_linked_distRL(opt1_events,opt2_events,choice,start,abandontn,resp_made,alphabins,biasbins,varsbins,fig_yes)
%inputs opt1_events[nblk,ntrial],opt2_events[nblk,ntrial],choice[nblk,ntrial]
if(nargin<10)fig_yes=0; end
if(nargin<9) varsbins=30; end
if(nargin<8) biasbins=30; end
if(nargin<7) alphabins=30; end
if(nargin<6) resp_made=true(size(choice)); end
if(nargin<5) abandontn=0; end
if(nargin<4) start=0.5; end

out=struct;

logit_alpha_space=logit(0.01):(logit(0.99)-logit(0.01))/(alphabins-1):logit(0.99);
logit_pos_bias_space=logit(0.02):(logit(0.98)-logit(0.02))/(biasbins-1):logit(0.98);
%log_vars_space=log(0.001):(log(0.3)-log(0.001))/(varsbins-1):log(0.3);
s_vars_space=0.001^2:(0.20^2-0.001^2)/(varsbins-1):0.20^2;

alpha_space=inv_logit(logit_alpha_space);
vars_space=sqrt(s_vars_space);
pos_bias_space=inv_logit(logit_pos_bias_space);
%vars_space=exp(s_vars_space);
for i=1:length(pos_bias_space)
    pos_bias=pos_bias_space(i);
    for j=1:length(vars_space)
        vars=vars_space(j);
        betab(j,i)=(pos_bias - vars + pos_bias.*vars - pos_bias.^2*2 + pos_bias.^3)./vars;
        betaa(j,i)=-(pos_bias.*(pos_bias.^2 - pos_bias + vars))./vars;
    end
end

for i=1:size(betaa,2)
    idx=find(betaa(:,i)<0);
    betaa(idx,i)=betaa(min(idx)-1,i);
    betab(idx,i)=betab(min(idx)-1,i);
end

for i=1:length(alpha_space)
    for j=1:length(vars_space)
        for k=1:length(pos_bias_space)
            [expv1,expv2]=cal_EV_distRL(alpha_space(i),betaa(j,k),betab(j,k),opt1_events,opt2_events,start);
            for m=1:size(expv1,2)    
                for l=abandontn+1:size(expv1,3)
                        probleft(i,j,k,m,l-abandontn)=read_prob(expv1(:,m,l),expv2(:,m,l));
                end
            end
        end
    end
end
choice=choice(:,abandontn+1:end);
ch=repmat(permute(choice,[5 3 4 1 2]),[length(alpha_space) length(vars_space) length(pos_bias_space) 1 1]);
% This calculates the likelihood of the choices made given the model
% parameters
probch=((ch.*probleft)+((1-ch).*(1-probleft)));
clear probleft

  %this bit removes data from trials in which the participant made no
  %response
 
  %probch=probch(:,:,resp_made,:,:);

  % This calculates the overall liklihood of the parameters by taking the
  % product of the individual trials. I presume the final term which
  % multiples the number by a large amount just makes the numbers
  % manageable (otherwise they are very small). Note that this is now a
  % four dimensional matrix which contains the likelihood of the data
  % given the three parameters which are coded on the dimensions (learning
  % rate, temperature, a). The final dimension separates the data from each
  % indvidiual participant
  out.posterior_prob(:,:,:)=squeeze(prod(prod(probch,4),5))*10^(size(probch,4)*size(probch,5)/5);  
  %renormalise
  out.posterior_prob=out.posterior_prob./(sum(sum(sum(out.posterior_prob))));
  
  clear probch


% This produces a marginal distribution of the learning rate by summing
% across the other dimensions and then renormalising.
% The numerator sums over dimensions 2 and 3 leaving a 2D matrix which contains
% the marginal likelihood for each of the different value of learning rate, for each
% participant. The denominator calcualtes the total likelihood for each of the subjects
% and then produces a matrix which reproduces these totals in each column
% with alphabins (length(i)) number of rows.
out.marg_base_alpha=squeeze(sum(sum(out.posterior_prob,3),2));

% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_base_alpha=inv_logit(logit_alpha_space*out.marg_base_alpha);

% this calculates the variance of the distribution of learning rates

out.var_base_alpha=inv_logit(((logit_alpha_space-logit(out.mean_base_alpha)).^2)*out.marg_base_alpha);


out.marg_vars=squeeze(sum(sum(out.posterior_prob,1),3))';
% out.mean_vars=exp(s_vars_space*out.marg_vars);
% out.var_vars=exp(((s_vars_space-log(out.mean_vars)).^2)*out.marg_vars);
out.mean_vars=sqrt(s_vars_space*out.marg_vars);
out.var_vars=sqrt(((s_vars_space-(out.mean_vars)).^4)*out.marg_vars);

out.marg_pos_bias=squeeze(sum(sum(out.posterior_prob,1),2));
out.mean_pos_bias=inv_logit(logit_pos_bias_space*out.marg_pos_bias);
out.var_pos_bias=inv_logit(((logit_pos_bias_space-logit(out.mean_pos_bias)).^2)*out.marg_pos_bias);

    

out.baselr_label=alpha_space;
out.vars_label=vars_space;
out.pos_bias_label=pos_bias_space;
out.mean_betab=(out.mean_pos_bias - out.mean_vars + out.mean_pos_bias.*out.mean_vars - out.mean_pos_bias.^2*2 + out.mean_pos_bias.^3)./out.mean_vars;
out.mean_betaa=-(out.mean_pos_bias.*(out.mean_pos_bias.^2 - out.mean_pos_bias + out.mean_vars))./out.mean_vars;

if fig_yes==1 
   figure
   subplot(3,1,1);
   plot(alpha_space,out.marg_base_alpha);
   title('base learning rate');
   subplot(3,1,2);
   plot(vars_space,out.marg_vars);
   title('variance');
   subplot(3,1,3);
   plot(pos_bias_space,out.marg_pos_bias);
    title('positve bias');

end

% get LL from estimated parameters
[expv1,expv2]=cal_EV_distRL(out.mean_base_alpha,out.mean_betaa,out.mean_betab,opt1_events,opt2_events,start,100);
for m=1:size(expv1,2)
    for l=abandontn+1:size(expv1,3)
        prob_ch_left(m,l-abandontn)=read_prob(expv1(:,m,l),expv2(:,m,l));
    end
end
resp_made=resp_made(:,abandontn+1:end);

likelihood=prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
out.neg_log_like=-sum(log(likelihood(resp_made)+1e-16));
out.BIC=(2.*out.neg_log_like)+3*(log(sum(sum(resp_made)))-log(2*pi)); % BIC given 3 free parameters
out.AIC=(2.*out.neg_log_like)+6;