function [out]=fit_linked_Bayesian_CVaR_eta_priorvar(opt1_events,opt2_events,choice,abandontn,resp_made,alpha_bins,eta_bins,var_bins,fig_yes)
%inputs opt1_events[nblk,ntrial],opt2_events[nblk,ntrial],choice[nblk,ntrial]
if(nargin<9) fig_yes=0; end
if(nargin<8) var_bins=30; end
if(nargin<7) eta_bins=30; end
if(nargin<6) alpha_bins=30; end
if(nargin<5) resp_made=true(size(choice)); end
if(nargin<4) abandontn=0; end

    out=struct;
    
    %eta_space=-0.95:(0.95+0.95)/(etabins-1):0.95;
%     a=1:98/29:99;
%     for i=1:30;b(i)=prctile(randn(1000000,1),a(i));end
    eta_space=-0.99:0.99*2/(eta_bins-1):0.99;
    %eta_space=0.01:0.98/(eta_bins-1):0.99;
    %logit_eta_space=logit(eta_space);
    
    
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
    prior_mean=0.5;

    log_var_space=log(0.005):(log(0.08)-log(0.005))/(var_bins-1):log(0.08);
    var_space=exp(log_var_space);

    for j=1:var_bins
        priorb(j)=(prior_mean - var_space(j) + prior_mean *var_space(j) - prior_mean ^2*2 + prior_mean ^3)/var_space(j);
        priora(j)=-(prior_mean *(prior_mean ^2 - prior_mean  + var_space(j)))/var_space(j);
    end
%     figure;for i=1:30;fplot(@(x) betapdf(x,priora(i),priorb(i)),[0 1]);hold on;end
    %%
    %%
   N=1000;
   step=1.01;
   x=0.01:0.98/(N-1):0.99;

    for a=1:alpha_bins
        alpha=alpha_space(a);  
        for v=1:var_bins
        %    updatevar=kernel_space(k);
            for n=1:size(opt1_events,1)
                rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora(v),priorb(v),alpha,updatevar,N);
                rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora(v),priorb(v),alpha,updatevar,N);      
                for j=1:length(eta_space)             
                   eta=eta_space(j)*2-1;
                   [r(n,:,1),mse(n,:,1)]=readout_CVaR_std(squeeze(rrdist(n,:,:,1)),eta,step);
                   [r(n,:,2),mse(n,:,2)]=readout_CVaR_std(squeeze(rrdist(n,:,:,2)),eta,step);
                   for t=1:size(opt1_events,2)
                       prob1(n,t,j,a,v)=normcdf(r(n,t,1)-r(n,t,2),0,sqrt(mse(n,t,1)+mse(n,t,2)));
                   end
                end               
            end
        end
    end




    ch=repmat(permute(choice,[1 2 3 4 5]),[1 1 length(eta_space) length(alpha_space) length(var_space)]);
    % This calculates the likelihood of the choices made given the model
    % parameters
    probch=((ch.*prob1)+((1-ch).*(1-prob1)));
    clear prob1 rrdist r mse

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
out.var_eta=((logit_eta_space-logit((out.mean_eta+1)/2)).^2)*out.marg_eta*2;

if alpha_bins==1
    out.marg_alpha=1;
    out.mean_alpha=1;
    out.var_alpha=0;  
    
    out.marg_var=squeeze(sum(out.posterior_prob,1))';
    out.mean_var=exp(log_var_space*out.marg_var);
    out.var_var=exp(((log_var_space-log(out.mean_var)).^2)*out.marg_var);
    
else
    out.marg_alpha=squeeze(sum(sum(out.posterior_prob,1),3))';
    out.mean_alpha=alpha_space*out.marg_alpha;
    out.var_alpha=(alpha_space-out.mean_alpha).^2*out.marg_alpha;
    
    out.marg_var=squeeze(sum(sum(out.posterior_prob,1),2));
    out.mean_var=exp(log_var_space*out.marg_var);
    out.var_var=exp(((log_var_space-log(out.mean_var)).^2)*out.marg_var);
    
end
  


out.eta_label=eta_space;
out.alpha_label=alpha_space;
out.var_label=var_space;

if fig_yes==1 
   figure
   subplot(1,3,1);
   plot(alpha_space,out.marg_alpha);
   title('alpha');
   subplot(1,3,2);
   plot(var_space,out.marg_var);
   title('var');
   subplot(1,3,3);
   plot(eta_space,out.marg_eta);
   title('eta');
end

% get LL from estimated parameters
priorb=(prior_mean - out.mean_var + prior_mean *out.mean_var - prior_mean ^2*2 + prior_mean ^3)/out.mean_var;
priora=-(prior_mean *(prior_mean ^2 - prior_mean  + out.mean_var))/out.mean_var;
for n=1:size(opt1_events,1)
    rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);
    rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);    
    [r(n,:,1),mse(n,:,1)]=readout_CVaR_std(squeeze(rrdist(n,:,:,1)),out.mean_eta,step);
    [r(n,:,2),mse(n,:,2)]=readout_CVaR_std(squeeze(rrdist(n,:,:,2)),out.mean_eta,step);
    for t=1:size(opt1_events,2)
        prob_ch_left(n,t)=normcdf(r(n,t,1)-r(n,t,2),0,sqrt(mse(n,t,1)+mse(n,t,2)));
    end
end


    


resp_made=resp_made(:,abandontn+1:end);

likelihood=prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
out.neg_log_like=-sum(log(likelihood(resp_made)+1e-16));
if alpha_bins==1
    out.BIC=(2.*out.neg_log_like)+2*log(sum(sum(resp_made))); % -2ln(LL)+kln(n)
    out.AIC=(2.*out.neg_log_like)+4;
else
    out.BIC=(2.*out.neg_log_like)+3*log(sum(sum(resp_made))); % -2ln(LL)+kln(n)
    out.AIC=(2.*out.neg_log_like)+4;
end