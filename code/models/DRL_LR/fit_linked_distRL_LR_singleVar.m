function [out]=fit_linked_distRL_LR_singleVar(opt1_events,opt2_events,choice,start,abandontn,resp_made,alpha_pos_bins,alpha_neg_bins,varsbins,fig_yes)
%inputs opt1_events[nblk,ntrial],opt2_events[nblk,ntrial],choice[nblk,ntrial]
if(nargin<10)fig_yes=0; end
if(nargin<9) varsbins=30; end
if(nargin<8) alpha_neg_bins=30; end
if(nargin<7) alpha_pos_bins=30; end
if(nargin<6) resp_made=true(size(choice)); end
if(nargin<5) abandontn=0; end
if(nargin<4) start=0.5; end

out=struct;

logit_alpha_pos_space=logit(0.2):(logit(0.8)-logit(0.2))/(alpha_pos_bins-1):logit(0.8);
logit_alpha_neg_space=logit(0.2):(logit(0.8)-logit(0.2))/(alpha_neg_bins-1):logit(0.8);
s_vars_space=0.001^2:(0.25^2-0.001^2)/(varsbins-1):0.25^2;

alpha_pos_space=inv_logit(logit_alpha_pos_space);
alpha_neg_space=inv_logit(logit_alpha_neg_space);
vars_space=sqrt(s_vars_space);

for i=1:length(alpha_pos_space)
    for j=1:length(alpha_neg_space)
        for k=1:length(vars_space)
            [expv1,expv2]=cal_EV_distRL_LR(alpha_pos_space(i),alpha_neg_space(j),vars_space(k),vars_space(k),opt1_events,opt2_events,start);
            for m=1:size(expv1,2)    
                for l=abandontn+1:size(expv1,3)
                        probleft(i,j,k,m,l-abandontn)=read_prob(expv1(:,m,l),expv2(:,m,l));
                end
            end
        end
    end
end
choice=choice(:,abandontn+1:end);
ch=repmat(permute(choice,[5 3 4 1 2]),[length(alpha_pos_space) length(alpha_neg_space) length(vars_space) 1 1]);
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
out.marg_pos_alpha=squeeze(sum(sum(out.posterior_prob,3),2));

% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_pos_alpha=inv_logit(logit_alpha_pos_space*out.marg_pos_alpha);

% this calculates the variance of the distribution of learning rates

out.var_pos_alpha=inv_logit(((logit_alpha_pos_space-logit(out.mean_pos_alpha)).^2)*out.marg_pos_alpha);


out.marg_neg_alpha=squeeze(sum(sum(out.posterior_prob,1),3))';
out.mean_neg_alpha=inv_logit(logit_alpha_neg_space*out.marg_neg_alpha);
out.var_neg_alpha=inv_logit(((logit_alpha_neg_space-logit(out.mean_neg_alpha)).^2)*out.marg_neg_alpha);

out.marg_vars=squeeze(sum(sum(out.posterior_prob,1),2));
out.mean_vars=sqrt(s_vars_space*out.marg_vars);
out.var_vars=sqrt(((s_vars_space-(out.mean_vars)).^4)*out.marg_vars);

    

out.alpha_pos_label=alpha_pos_space;
out.vars_label=vars_space;
out.alpha_neg_label=alpha_neg_space;
out.betab_plr=(out.mean_pos_alpha - out.mean_vars + out.mean_pos_alpha.*out.mean_vars - out.mean_pos_alpha.^2*2 + out.mean_pos_alpha.^3)./out.mean_vars;
out.betaa_plr=-(out.mean_pos_alpha.*(out.mean_pos_alpha.^2 - out.mean_pos_alpha + out.mean_vars))./out.mean_vars;

out.betab_nlr=(out.mean_neg_alpha - out.mean_vars + out.mean_neg_alpha.*out.mean_vars - out.mean_neg_alpha.^2*2 + out.mean_neg_alpha.^3)./out.mean_vars;
out.betaa_nlr=-(out.mean_neg_alpha.*(out.mean_neg_alpha.^2 - out.mean_neg_alpha + out.mean_vars))./out.mean_vars;

out.pos_bias=out.mean_pos_alpha/(out.mean_neg_alpha+out.mean_pos_alpha);
out.baselrs=(out.mean_neg_alpha+out.mean_pos_alpha)/2;
if fig_yes==1 
   % for i=1:94;
    %    out=estimates_DRL_LR(i);
   figure
   subplot(3,1,1);
   plot(alpha_pos_space,out.marg_pos_alpha);
   title('positive learning rate');
   subplot(3,1,2);
   plot(vars_space,out.marg_vars);
   title('variance');
   subplot(3,1,3);
   plot(alpha_neg_space,out.marg_neg_alpha);
    title('negative learning rate');
   % end

end

% get LL from estimated parameters
[expv1,expv2]=cal_EV_distRL_LR(out.mean_pos_alpha,out.mean_neg_alpha,out.mean_vars,out.mean_vars,opt1_events,opt2_events,start,100);
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