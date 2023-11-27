function [out]=fit_linked_Bayesian_CVaR_eta_beta(opt1_events,opt2_events,choice,abandontn,resp_made,alpha_bins,eta_bins,betabins,fig_yes)
%inputs opt1_events[nblk,ntrial],opt2_events[nblk,ntrial],choice[nblk,ntrial]
if(nargin<9) fig_yes=0; end
if(nargin<8) betabins=30; end
if(nargin<7) eta_bins=30; end
if(nargin<6) alpha_bins=30; end
if(nargin<5) resp_made=true(size(choice)); end
if(nargin<4) abandontn=0; end

    out=struct;
    %log_beta_space=log(0.01):(log(20)-log(0.01))/(betabins-1):log(20);
    log_beta_space=log(0.05):(log(20)-log(0.05))/(betabins-1):log(20);
    beta_space=exp(log_beta_space);
%    figure;
%     for i=1:length(beta_space)
%         fplot(@(x) 1./(1+exp(-beta_space(i)*(x))),[-1 1]);
%         hold on;
%     end
    %eta_space=-0.95:(0.95+0.95)/(etabins-1):0.95;
%     a=1:98/29:99;
%     for i=1:30;b(i)=prctile(randn(1000000,1),a(i));end
%     tmp=betarnd(2,2,(eta_bins-1)*1000+2*eta_bins,1);
%     tmp=sort(tmp);
%     eta_space=transpose(tmp(eta_bins:1000:(eta_bins-1)*1000+eta_bins));
    
%     eta_space=0.01:0.98/(eta_bins-1):0.99;
%     logit_eta_space=logit(eta_space);
    logit_eta_space=logit(0.001):(logit(0.999)-logit(0.001))/(eta_bins-1):logit(0.999);   
    eta_space=inv_logit(logit_eta_space);

    if alpha_bins==1
        alpha_space=1;
    else
        logit_alpha_space=logit(0.01):(logit(0.99)-logit(0.01))/(alpha_bins-1):logit(0.99);
        alpha_space=inv_logit(logit_alpha_space);
    end
    %kernel_bins=20;
    updatevar=0.009;%updatevar=0.009; %updatevar=0.0005;%
    % log_kernel_space=log(0.0001):(log(0.0099)-log(0.0001))/(kernel_bins-1):log(0.0099);
    % kernel_space=exp(log_kernel_space);
    %
%     prior_mean=0.5;
%     prior_var=0.06;
%     priorb=(prior_mean - prior_var + prior_mean *prior_var - prior_mean ^2*2 + prior_mean ^3)/prior_var;
%     priora=-(prior_mean *(prior_mean ^2 - prior_mean  + prior_var))/prior_var;

    priorb=1;
    priora=1;
%     figure;fplot(@(x) betapdf(x,priora,priorb),[0 1])
    %%
    %%
   N=10000;
   step=1.01;


    for a=1:alpha_bins
        alpha=alpha_space(a);
        %for k=1:kernel_bins
        %    updatevar=kernel_space(k);
            for n=1:size(opt1_events,1)
                rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora,priorb,alpha,updatevar,N);
                rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora,priorb,alpha,updatevar,N);           
                for j=1:length(eta_space)             
                   eta=eta_space(j)*2-1;
                   r(n,:,1)=readout_CVaR(squeeze(rrdist(n,:,:,1)),eta,step);
                   r(n,:,2)=readout_CVaR(squeeze(rrdist(n,:,:,2)),eta,step);
                   for i=1:betabins
                       beta=beta_space(i);
                       prob1(n,:,j,a,i)=squeeze(exp(beta*r(n,:,1))./(exp(beta*r(n,:,1))+exp(beta*r(n,:,2))));
                   end
                end
            end
        %end
    end





    ch=repmat(permute(choice,[1 2 3 4 5]),[1 1 length(eta_space) length(alpha_space) length(beta_space)]);
    % This calculates the likelihood of the choices made given the model
    % parameters
    probch=((ch.*prob1)+((1-ch).*(1-prob1)));
    clear prob1 rrdist r

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
  out.posterior_prob(:,:,:)=squeeze(prod(prod(probch,1),2))*10^(size(probch,1)*size(probch,2)/6);  
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


% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
% this calculates the variance of the distribution of learning rates
out.marg_eta=squeeze(sum(sum(out.posterior_prob,2),3));
out.mean_eta=inv_logit(logit_eta_space*out.marg_eta)*2-1;
%out.mean_eta=eta_space*out.marg_eta;
out.var_eta=((logit_eta_space-logit((out.mean_eta+1)/2)).^2)*out.marg_eta*2;
%out.var_eta=(eta_space-out.mean_eta).^2*out.marg_eta;

if alpha_bins==1
    out.marg_alpha=1;
    out.mean_alpha=1;
    out.var_alpha=0;  
    
    out.marg_beta=squeeze(sum(out.posterior_prob,1))';
    out.mean_beta=exp(log_beta_space*out.marg_beta);
    out.var_beta=exp(((log_beta_space-log(out.mean_beta)).^2)*out.marg_beta);
    
else
    out.marg_alpha=squeeze(sum(sum(out.posterior_prob,1),3))';
    out.mean_alpha=alpha_space*out.marg_alpha;
    out.var_alpha=(alpha_space-out.mean_alpha).^2*out.marg_alpha;
    
    out.marg_beta=squeeze(sum(sum(out.posterior_prob,1),2));
    out.mean_beta=exp(log_beta_space*out.marg_beta);
    out.var_beta=exp(((log_beta_space-log(out.mean_beta)).^2)*out.marg_beta);
    
end
  


out.eta_label=eta_space*2-1;
out.alpha_label=alpha_space;
out.beta_label=beta_space;

if fig_yes==1 
   figure
   subplot(1,3,1);
   plot(alpha_space,out.marg_alpha);
   title('alpha');
   subplot(1,3,2);
   plot(beta_space,out.marg_beta);
   title('beta');
   subplot(1,3,3);
   plot(eta_space,out.marg_eta);
   title('eta');
end

% get LL from estimated parameters
for n=1:size(opt1_events,1)
    rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);
    rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);    
    r(n,:,1)=readout_CVaR(squeeze(rrdist(n,:,:,1)),out.mean_eta,step);
    r(n,:,2)=readout_CVaR(squeeze(rrdist(n,:,:,2)),out.mean_eta,step);
end
    %for j=1:length(eta_space)             
    %  eta=eta_space(j);
    out.prob_ch_left=squeeze(exp(out.mean_beta*r(:,:,1))./(exp(out.mean_beta*r(:,:,1))+exp(out.mean_beta*r(:,:,2))));
    


resp_made=resp_made(:,abandontn+1:end);

likelihood=out.prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
out.neg_log_like=-sum(log(likelihood(resp_made)+1e-16));
if alpha_bins==1
    out.BIC=(2.*out.neg_log_like)+2*log(sum(sum(resp_made))); % -2ln(LL)+kln(n)
    out.AIC=(2.*out.neg_log_like)+4;
else
    out.BIC=(2.*out.neg_log_like)+3*log(sum(sum(resp_made))); % -2ln(LL)+kln(n)
    out.AIC=(2.*out.neg_log_like)+4;
end