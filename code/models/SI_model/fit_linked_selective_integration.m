function [out]=fit_linked_selective_integration(opt1_events,opt2_events,choice,start,abandontn,resp_made,wbins,leakbins,noisebins,fig_yes)

if(nargin<10)fig_yes=0; end
if(nargin<9) noisebins=30; end
if(nargin<8) leakbins=30; end
if(nargin<7) wbins=30; end
if(nargin<6) resp_made=true(size(choice)); end
if(nargin<5) abandontn=0; end
if(nargin<4) start=0.5; end

out=struct;
maxb=3;
maxb1=3;
log_w_space=log(0.01):(log(maxb)-log(0.01))/(wbins-1):log(maxb);
w_space=exp(log_w_space);
leak_space=0:0.5/(leakbins-1):0.5;
log_noise_space=log(0.01):(log(maxb1)-log(0.01))/(noisebins-1):log(maxb1);
noise_space=exp(log_noise_space);
lapse_space=0:0.2/(noisebins-1):0.2;

for i=1:length(w_space)
    for j=1:length(leak_space)
        for k=1:length(noise_space)
            for l=1:length(lapse_space)
                 probleft(i,j,k,l,:,:)=selective_gating(opt1_events,opt2_events,w_space(i),leak_space(j),noise_space(k),lapse_space(l),abandontn,start);       
            end
        end
    end
end

choice=choice(:,abandontn+1:end);  
% create a 5D matrix of the same dimensions as before with the choices
% arranged along the second dimension
ch=repmat(permute(choice,[6 4 5 3 1 2]),[length(w_space) length(leak_space) length(noise_space) length(lapse_space) 1 1]);
% This calculates the likelihood of the choices made given the model
% parameters
probch=((ch.*probleft)+((1-ch).*(1-probleft)));
clear probleft

  %this bit removes data from trials in which the participant made no
  %response
 
  %probch=probch(:,:,resp_made,:,:);
 % includeSample=[false(abandontn,size(opt1_events,2));resp_made((abandontn+1):end,:)];
 
%   for i=1:size(opt1_events,2)
%   tmp=probch(:,:,:,:,includeSample(:,1),i);
  out.posterior_prob(:,:,:,:)=squeeze(prod(prod(probch,5),6))*10^(size(probch,5)*size(probch,6)/5);  
 % end 
  
%  clear tmp
   clear probch
%   out.posterior_prob=squeeze(prod(posterior_probt,5))*10^(size(posterior_probt,5)/5); 
%   clear posterior_probt
  % This calculates the overall liklihood of the parameters by taking the
  % product of the individual trials. I presume the final term which
  % multiples the number by a large amount just makes the numbers
  % manageable (otherwise they are very small). Note that this is now a
  % four dimensional matrix which contains the likelihood of the data
  % given the three parameters which are coded on the dimensions (learning
  % rate, temperature, a). The final dimension separates the data from each
  % indvidiual participant
  %renormalise
  out.posterior_prob=out.posterior_prob./(sum(sum(sum(sum(out.posterior_prob)))));


% This produces a marginal distribution of the learning rate by summing
% across the other dimensions and then renormalising.
% The numerator sums over dimensions 2 and 3 leaving a 2D matrix which contains
% the marginal likelihood for each of the different value of learning rate, for each
% participant. The denominator calcualtes the total likelihood for each of the subjects
% and then produces a matrix which reproduces these totals in each column
% with wbins (length(i)) number of rows.
out.marg_w=squeeze(sum(sum(sum(out.posterior_prob,3),2),4));

% This generates the logected value of the learning rate using a weighted
% sum-- marginal probabilities multiplied by learning rate values. Note
% for both the learning rate and temperature mean and variance are
% caculated in log space
out.mean_w=exp(log_w_space*out.marg_w);

% this calculates the variance of the distribution of learning rates

out.var_w=exp((log_w_space-log(out.mean_w)).^2*out.marg_w);


out.marg_leak=squeeze(sum(sum(sum(out.posterior_prob,1),3),4))';
out.mean_leak=leak_space*out.marg_leak;
out.var_leak=(leak_space-out.mean_leak).^2*out.marg_leak;

out.marg_noise=squeeze(sum(sum(sum(out.posterior_prob,1),2),4));
out.mean_noise=exp(log_noise_space*out.marg_noise);
out.var_noise=exp((log_noise_space-log(out.mean_noise)).^2*out.marg_noise);

out.marg_lapse=squeeze(sum(sum(sum(out.posterior_prob,1),2),3));
out.mean_lapse=lapse_space*out.marg_lapse;
out.var_lapse=(lapse_space-out.mean_lapse).^2*out.marg_lapse;
    
out.noise_label=noise_space;
out.w_label=w_space;
out.leak_label=leak_space;
out.lapse_space=lapse_space;

if fig_yes==1 
   figure
   subplot(2,2,1);
   plot(w_space,out.marg_w);
   title('w');
   subplot(2,2,2);
   plot(leak_space,out.marg_leak);
   title('leak');
   subplot(2,2,3);
   plot(noise_space,out.marg_noise);
    title('early noise');
      subplot(2,2,4);
    plot(lapse_space,out.marg_lapse);
    title('lapse');   
end

% get LL from estimated parameters
prob_ch_left=selective_gating(opt1_events,opt2_events,out.mean_w,out.mean_leak,out.mean_noise,out.mean_lapse,abandontn,start);
likelihood=prob_ch_left;
likelihood(choice==0)=1-likelihood(choice==0);
resp_made=resp_made(:,abandontn+1:end);
out.neg_log_like=-sum(log(likelihood(resp_made)+1e-16));
out.BIC=(2.*out.neg_log_like)+4*(log(sum(sum(resp_made)))-log(2*pi)); % BIC given 4 free parameters
out.AIC=(2.*out.neg_log_like)+8;