clear all
clc
load('simulated_data_SI_model.mat')
abandontn=5;
parfor i=1:size(choices,1)%48
    disp(i)
    choice=squeeze(choices(i,:,:));
    estimates=fit_linked_selective_integration(squeeze(blocks_events(:,1,:)),squeeze(blocks_events(:,2,:)),choice,abandontn);
    pararecovered(i,:)=[estimates.mean_w,estimates.mean_leak,estimates.mean_noise,estimates.mean_lapse];
end
save('recovered_data_SI.mat')

variablenames={'w','leak','noise','lapse'};

 for i=1:2:4
 f=figure;scatter(log(paraset(1:100,i)),log(pararecovered(:,i)))
 lsline
 title(variablenames(i))
 xlabel('true values -log')
 ylabel('recovered values -log')
 saveas(f,['recovery_',variablenames{i},'.png'])
 end
  for i=2:2:5
 f=figure;scatter(paraset(1:100,i),pararecovered(:,i))
 lsline
 title(variablenames(i))
 xlabel('true values')
 ylabel('recovered values')
 saveas(f,['recovery_',variablenames{i},'.png'])
 end
% estimates=fit_linked_selective_integration(opt1_events(trialrange),opt2_events(trialrange),choice,abandontn,resp_made,wbins,leakbins,noisebins,1);