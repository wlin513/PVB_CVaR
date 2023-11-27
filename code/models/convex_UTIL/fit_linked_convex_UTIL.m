function [out]=fit_linked_convex_UTIL(opt1_events,opt2_events,choice,start,abandontn,resp_made,alphabins,betabins,kbins,fig_yes)


% Bayesian model fit for 2 option vol training task in which outcomes are linked.
% (i.e. if one option wins, the other doesn't).
% Information contains left_win, left_loss (values for right are 1-these)
% choice is the participant choice, start [rew_start loss_start]
% alphabins sets number of points to
% estimate lr, betabins same for beta, resp_made is trials on which a
% response was made.
% in this version different betas independently multiply win and loss values 
if(nargin<10) fig_yes=0; end
if(nargin<9) kbins=30; end 
if(nargin<8) betabins=30; end 
if(nargin<7) alphabins=30; end 
if(nargin<6) resp_made=(isnan(choice)==0); end  %exclude the trials on which no response is made

out=struct;

% NB calculations (of mean and variance) of both learning rate and decision
% temperature are performed in log space.

%This creates a vector of length alphabins in inv_logit space
i=logit(0.01):(logit(0.99)-logit(0.01))/(alphabins-1):logit(0.99);

% this creates a vector of length betabins the value of which changes
% linearly from log(1)=0 to log(100)=4.6. It will be used to create a
% logorythmic distribution of inverse temperatures.
b_label=log(0.05):(log(20)-log(0.05))/(betabins-1):log(20);

%k_space=-2:4/(kbins-1):2;
 k_space=log(0.2):(log(8)-log(0.2))/(kbins-1):log(8);
% figure;
% for kk=1:20%length(k_space)
%     fplot(@(x) exponential_utility(x,-exp(k_space(kk))),[0 1]);
%     hold on
% end
%   fplot(@(x) x,[0 1]);

for a=1:length(i)
    for kk=1:length(k_space)
        learn_left=rescorla_wagner_UTIL(opt1_events,opt2_events,inv_logit(i(a)),-exp(k_space(kk)),start);
        val_l(a,kk,:)=learn_left; 
    end
end

% replicate the value matrice to account for betas

mmdl=repmat(val_l,[1,1,1,length(b_label)]);

clear val_l

beta=repmat(permute(exp(b_label),[1 3 4 2]),[length(i) length(k_space) size(opt1_events,1)*size(opt2_events,2) 1]);

probleft=1./(1+exp(-(beta.*mmdl)));

clear beta
% looks like choices can be a matrix with each column representing the
% choice of a given participant


  
  % create a 5D matrix of the same dimensions as before with the choices
  % arranged along the second dimension
  chtmp=reshape(choice',size(choice,1)*size(choice,2),1);
  ch=repmat(permute(chtmp,[3 4 1 2]),[length(i) length(k_space) 1 length(b_label)]);
  % This calculates the likelihood of the choices made given the model
  % parameters
  probch=((ch.*probleft)+((1-ch).*(1-probleft)));
  clear probleft

  %this bit removes data from trials in which the participant made no
  %response
  includeTri=resp_made;
  includeTri(:,1:abandontn)=false(size(opt1_events,1),abandontn);
  includetrials=reshape(includeTri',size(includeTri,1)*size(includeTri,2),1);
  probch=probch(:,:,includetrials,:);


  % This calculates the overall liklihood of the parameters by taking the
  % product of the individual trials. I presume the final term which
  % multiples the number by a large amount just makes the numbers
  % manageable (otherwise they are very small). Note that this is now a
  % four dimensional matrix which contains the likelihood of the data
  % given the three parameters which are coded on the dimensions (learning
  % rate, temperature, a). The final dimension separates the data from each
  % indvidiual participant
  out.posterior_prob(:,:,:)=squeeze(prod(probch,3))*10^(size(probch,2)/5);  
  %renormalise
  out.posterior_prob=out.posterior_prob./(sum(sum(sum(out.posterior_prob))));
  
  clear probch

% this returns the actual values of the parameters used (I guess for
% graphing)
alphalabel=inv_logit(i);
betalabel=exp(b_label);
klabel=exp(k_space);

% This produces a marginal distribution of the learning rate by summing
% across the other dimensions and then renormalising.
% The numerator sums over dimensions 2 and 3 leaving a 2D matrix which contains
% the marginal likelihood for each of the different value of learning rate, for each
% participant. The denominator calcualtes the total likelihood for each of the subjects
% and then produces a matrix which reproduces these totals in each column
% with alphabins (length(i)) number of rows.
out.marg_alpha=squeeze(sum(sum(out.posterior_prob,2),3));

% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_alpha=inv_logit(i*out.marg_alpha);

% this calculates the variance of the distribution of learning rates

out.var_alpha=inv_logit(((i-logit(out.mean_alpha)).^2)*out.marg_alpha);

    out.marg_k=squeeze(sum(sum(out.posterior_prob,1),3))';
    out.mean_k=-exp(k_space*out.marg_k);
    out.var_k=exp(((k_space-log(out.mean_k)).^2)*out.marg_k);

%beta 

    out.marg_beta=squeeze(sum(sum(out.posterior_prob,1),2));
    out.mean_beta=exp(b_label*out.marg_beta);
    out.var_beta=exp(((b_label-log(out.mean_beta)).^2)*out.marg_beta);
    

out.beta_label=betalabel;
out.lr_label=alphalabel;
out.lr_points=i;
out.beta_points=b_label;
out.k_label=klabel;



if fig_yes==1 
   figure
   subplot(3,1,1);
   plot(alphalabel,out.marg_alpha);
   title('Learning Rate');
   subplot(3,1,2);
   plot(betalabel,out.marg_beta);
    title('Beta');

    subplot(3,1,3);
   plot(klabel,out.marg_k);
    title('k');

   
end
% get LL from estimated parameters
bel=rescorla_wagner_UTIL(opt1_events,opt2_events,out.mean_alpha,out.mean_k,start);
out.prob_ch_left=1./(1+exp(-(out.mean_beta.*bel)));
likelihood=out.prob_ch_left;
likelihood(chtmp==0)=1-likelihood(chtmp==0);
out.neg_log_like=-sum(log(likelihood(includetrials)+1e-16));
out.BIC=(2.*out.neg_log_like)+3*log(sum(sum(resp_made))); % -2ln(LL)+kln(n)
out.AIC=(2.*out.neg_log_like)+6;
   