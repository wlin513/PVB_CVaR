function [out]=fit_linked_PEIRS_fix_s0(opt1_events,opt2_events,choice,startQ,startS,abandontn,resp_made,alphabins,betabins,omegabins,fig_yes)


% Bayesian model fit
if(nargin<11) fig_yes=0; end
if(nargin<10) betabins=30; end 
if(nargin<9) alphabins=30; end 
if(nargin<8) omegabins=30; end 
if(nargin<7) resp_made=(isnan(choice)==0); end  %exclude the trials on which no response is made

out=struct;

% NB calculations (of mean and variance) of both learning rate and decision
% temperature are performed in log space.

%This creates a vector of length alphabins in inv_logit space
i=logit(0.01):(logit(0.99)-logit(0.01))/(alphabins-1):logit(0.99);

% this creates a vector of length betabins the value of which changes
% linearly from log(1)=0 to log(100)=4.6. It will be used to create a
% logorythmic distribution of inverse temperatures.
b_log=log(0.1):(log(100)-log(0.1))/(betabins-1):log(100);

orange=5;
o=-orange:(orange-(-orange))/(omegabins-1):orange;

for k=1:length(i)
for j=1:length(i)
    [v1(k,j,:,:),v2(k,j,:,:),s1(k,j,:,:),s2(k,j,:,:)]=PEIRS(opt1_events,opt2_events,inv_logit(i(k)),inv_logit(i(j)),startQ,startS);

end
end
pes=(v1+v2)./2-0.5;%stimulus prediction errors
% replicate the value matrice to account for betas
for q=1:length(o)
    q1(:,:,:,:,q)=v1+tanh(o(q).*pes).*s1;
    q2(:,:,:,:,q)=v2+tanh(o(q).*pes).*s2;
end
clear v1 v2 s1 s2 pes
val_l=q1-q2;
mmdl=repmat(val_l,[1,1,1,1,1,length(b_log)]);

clear val_l q1 q2 

beta=repmat(permute(exp(b_log),[1 3 4 5 6 2]),[length(i) length(i) size(opt1_events,1) size(opt2_events,2) length(o) 1]);

probleft=1./(1+exp(-(beta.*mmdl)));

clear beta mmdl
% looks like choices can be a matrix with each column representing the
% choice of a given participant


  
  % create a 5D matrix of the same dimensions as before with the choices
  % arranged along the second dimension
  ch=repmat(permute(choice,[3 4 1 2 5 6]),[length(i) length(i) 1 1 length(o) length(b_log)]);
  % This calculates the likelihood of the choices made given the model
  % parameters
  probch=((ch.*probleft)+((1-ch).*(1-probleft)));
  clear probleft ch

  %this bit removes data from trials in which the participant made no
  %response

  probch(:,:,:,1:abandontn,:,:)=[];

  % This calculates the overall liklihood of the parameters by taking the
  % product of the individual trials. I presume the final term which
  % multiples the number by a large amount just makes the numbers
  % manageable (otherwise they are very small). Note that this is now a
  % four dimensional matrix which contains the likelihood of the data
  % given the three parameters which are coded on the dimensions (learning
  % rate, temperature, a). The final dimension separates the data from each
  % indvidiual participant
  out.posterior_prob=squeeze(prod(prod(probch,4),3))*10^(size(probch,4)*size(probch,3)/5);  
  %renormalise
  out.posterior_prob=out.posterior_prob./(sum(sum(sum(sum(out.posterior_prob)))));
  
  clear probch

% this returns the actual values of the parameters used (I guess for
% graphing)
alphalabel=inv_logit(i);
betalabel=exp(b_log);
omegalabel=o;

% This produces a marginal distribution of the learning rate by summing
% across the other dimensions and then renormalising.
% The numerator sums over dimensions 2 and 3 leaving a 2D matrix which contains
% the marginal likelihood for each of the different value of learning rate, for each
% participant. The denominator calcualtes the total likelihood for each of the subjects
% and then produces a matrix which reproduces these totals in each column
% with alphabins (length(i)) number of rows.
out.marg_alphaQ=squeeze(sum(sum(sum(out.posterior_prob,4),3),2));

% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_alphaQ=inv_logit(i*out.marg_alphaQ);

% this calculates the variance of the distribution of learning rates

out.var_alphaQ=inv_logit(((i-logit(out.mean_alphaQ)).^2)*out.marg_alphaQ);


out.marg_alphaS=squeeze(sum(sum(sum(out.posterior_prob,1),3),4))';

% This generates the expected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_alphaS=inv_logit(i*out.marg_alphaS);

% this calculates the variance of the distribution of learning rates

out.var_alphaS=inv_logit(((i-logit(out.mean_alphaS)).^2)*out.marg_alphaS);

    

out.marg_omega=squeeze(sum(sum(sum(out.posterior_prob,1),2),4));
out.mean_omega=o*out.marg_omega;
out.var_omega=(o-out.mean_omega).^2*out.marg_omega;
%beta 

    out.marg_beta=squeeze(sum(sum(sum(out.posterior_prob,1),2),3));
    out.mean_beta=exp(b_log*out.marg_beta);
    out.var_beta=exp(((b_log-log(out.mean_beta)).^2)*out.marg_beta);
    

out.beta_label=betalabel;
out.lr_label=alphalabel;
out.omegalabel=omegalabel;
out.lr_points=i;
out.beta_points=b_log;



if fig_yes==1 
   f=figure;
   subplot(2,2,1);
   plot(alphalabel,out.marg_alphaQ);
   title('Q Learning Rate');
   subplot(2,2,2);
   plot(alphalabel,out.marg_alphaS);
   title('S Learning Rate');
   subplot(2,2,3);
   plot(betalabel,out.marg_beta);
    title('Beta');

      subplot(2,2,4);
   plot(omegalabel,out.marg_omega);
    title('omega');
   
end
%out.figure=f;
% get LL from estimated parameters
[v1,v2,s1,s2]=PEIRS(opt1_events,opt2_events,out.mean_alphaQ, out.mean_alphaS,startQ,startS);
pes=(v1+v2)./2-0.5;
bel=v1+tanh(out.mean_omega.*pes).*s1-v2-tanh(out.mean_omega.*pes).*s2;
prob_ch_left=1./(1+exp(-(out.mean_beta.*bel)));
likelihood=prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
likelihood(:,1:abandontn)=[];
out.neg_log_like=-sum(sum(log(likelihood)+1e-16));
out.BIC=(2.*out.neg_log_like)+4*(log(size(likelihood,1)*size(likelihood,2))-log(2*pi)); % BIC given 4 free parameters
out.AIC=(2.*out.neg_log_like)+8;
   