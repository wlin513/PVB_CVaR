function [out]=fit_linked_Bayesian_CVaR(opt1_events,opt2_events,choice,abandontn,resp_made,alpha_bins,bias_bins,var_bins,betabins,fig_yes)
%inputs opt1_events[nblk,ntrial],opt2_events[nblk,ntrial],choice[nblk,ntrial]
if(nargin<10) fig_yes=0; end
if(nargin<9) betabins=30; end
if(nargin<8) var_bins=30; end
if(nargin<7) bias_bins=30; end
if(nargin<6) alpha_bins=30; end
if(nargin<5) resp_made=true(size(choice)); end
if(nargin<4) abandontn=0; end

    out=struct;

    log_beta_space=log(2):(log(50)-log(2))/(betabins-1):log(50);
    beta_space=exp(log_beta_space);

    %eta_space=-0.95:(0.95+0.95)/(etabins-1):0.95;

    logit_bias_space=logit(0.06):(logit(0.94)-logit(0.06))/(bias_bins-1):logit(0.94);
    bias_space=inv_logit(logit_bias_space);

    log_var_space=log(0.005):(log(0.05)-log(0.005))/(var_bins-1):log(0.05);
    var_space=exp(log_var_space);
    if alpha_bins==1
        alpha_space=1;
    else
        alpha_space=0:1/(alpha_bins-1):1;
    end
    
    updatevar=0.005;
%     log_kernel_space=log(0.0001):(log(0.0099)-log(0.0001))/(kernel_bins-1):log(0.0099);
%     kernel_space=exp(log_kernel_space);
    %
    for i=1:bias_bins
        for j=1:var_bins
            betab(i,j)=(bias_space(i) - var_space(j) + bias_space(i) *var_space(j) - bias_space(i) ^2*2 + bias_space(i) ^3)/var_space(j);
            betaa(i,j)=-(bias_space(i) *(bias_space(i) ^2 - bias_space(i)  + var_space(j)))/var_space(j);
        end
    end
%     figure;imagesc(betaa>0)
%     figure;for i=1:30;fplot(@(x) betapdf(x,betaa(30,i),betab(30,i)),[0 1]);hold on;end
    %%
    %%
   N=100;
   step=1.01;
   eta=0;
    for b=1:bias_bins
        for v=1:var_bins
            priora=betaa(b,v);
            priorb=betab(b,v);
            for a=1:alpha_bins
                alpha=alpha_space(a);
                %for k=1:kernel_bins
                    %updatevar=kernel_space(k);
                    for n=1:size(opt1_events,1)
                        rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora,priorb,alpha,updatevar,N);
                        rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora,priorb,alpha,updatevar,N);           
                        %for j=1:length(eta_space)             
                         %  eta=eta_space(j);
                           r(n,:,1)=readout_CVaR(squeeze(rrdist(n,:,:,1)),eta,step);
                           r(n,:,2)=readout_CVaR(squeeze(rrdist(n,:,:,2)),eta,step);
                           for i=1:betabins
                               beta=beta_space(i);
                               prob1(n,:,b,v,a,i)=squeeze(exp(beta*r(n,:,1))./(exp(beta*r(n,:,1))+exp(beta*r(n,:,2))));
                           end
                        %end
                    end
                %end
            end
        end
    end



    ch=repmat(permute(choice,[1 2 3 4 5 6]),[1 1 length(bias_space) length(var_space) length(alpha_space) length(beta_space)]);
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
  out.posterior_prob(:,:,:,:)=squeeze(prod(prod(probch,1),2))*10^(size(probch,1)*size(probch,2)/6);  
  %renormalise
  out.posterior_prob=out.posterior_prob./(sum(sum(sum(sum(out.posterior_prob)))));
  
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
out.marg_bias=squeeze(sum(sum(sum(out.posterior_prob,2),3),4));
out.mean_bias=inv_logit(logit_bias_space*out.marg_bias);
out.var_bias=inv_logit(logit_bias_space-logit(out.mean_bias).^2*out.marg_bias);

out.marg_var=squeeze(sum(sum(sum(out.posterior_prob,1),3),4))';
out.mean_var=exp(log_var_space*out.marg_var);
out.var_var=exp(((log_var_space-log(out.mean_var)).^2)*out.marg_var);

if alpha_bins==1
    out.marg_alpha=1;
    out.mean_alpha=1;
    out.var_alpha=0;  
    
    out.marg_beta=squeeze(sum(sum(out.posterior_prob,1),2));
    out.mean_beta=exp(log_beta_space*out.marg_beta);
    out.var_beta=exp(((log_beta_space-log(out.mean_beta)).^2)*out.marg_beta);
else
    out.marg_alpha=squeeze(sum(sum(sum(out.posterior_prob,1),2),4));
    out.mean_alpha=alpha_space*out.marg_alpha;
    out.var_alpha=(alpha_space-out.mean_alpha).^2*out.marg_alpha;
    
    out.marg_beta=squeeze(sum(sum(sum(out.posterior_prob,1),2),3));
    out.mean_beta=exp(log_beta_space*out.marg_beta);
    out.var_beta=exp(((log_beta_space-log(out.mean_beta)).^2)*out.marg_beta);
end



out.bias_label=bias_space;
out.var_label=var_space;
out.alpha_label=alpha_space;
out.beta_label=beta_space;

if fig_yes==1 
   figure
   subplot(2,2,1);
   plot(alpha_space,out.marg_alpha);
   title('alpha');
   subplot(2,2,2);
   plot(beta_space,out.marg_beta);
   title('beta');
   subplot(2,2,3);
   plot(bias_space,out.marg_bias);
   title('bias');
   subplot(2,2,4);
   plot(var_space,out.marg_var);
   title('var');
end

% get LL from estimated parameters
priorb=(out.mean_bias - out.mean_var + out.mean_bias *out.mean_var - out.mean_bias ^2*2 + out.mean_bias ^3)/out.mean_var;
priora=-(out.mean_bias *(out.mean_bias ^2 - out.mean_bias  + out.mean_var))/out.mean_var;
for n=1:size(opt1_events,1)
    rrdist(n,:,:,1)=bayesian_update(opt1_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);
    rrdist(n,:,:,2)=bayesian_update(opt2_events(n,:),priora,priorb,out.mean_alpha,updatevar,N);    
    r(n,:,1)=readout_CVaR(squeeze(rrdist(n,:,:,1)),eta,step);
    r(n,:,2)=readout_CVaR(squeeze(rrdist(n,:,:,2)),eta,step);
end
    %for j=1:length(eta_space)             
    %  eta=eta_space(j);
    prob_ch_left=squeeze(exp(out.mean_beta*r(:,:,1))./(exp(out.mean_beta*r(:,:,1))+exp(out.mean_beta*r(:,:,2))));
    


resp_made=resp_made(:,abandontn+1:end);

likelihood=prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
out.neg_log_like=-sum(log(likelihood(resp_made)+1e-16));
out.BIC=(2.*out.neg_log_like)+4*(log(sum(sum(resp_made)))-log(2*pi)); % BIC given 4 free parameters
out.AIC=(2.*out.neg_log_like)+8;