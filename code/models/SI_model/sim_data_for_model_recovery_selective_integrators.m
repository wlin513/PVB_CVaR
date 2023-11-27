clear
clc
%%  set up event distributions
load('../events.mat')
%1=Broad_high_event(BH_events);2=narrow_high_event(NH_events);
%3=Broad_low_event(BL_events);4=narrow_low_event(NL_events)
blocks_list=[1,2;3,4;1,4;1,3;2,3;2,4];
%blocks_list=[1,2;1,2;1,2;3,4;3,4;3,4;1,4;1,4;1,3;1,3;3,2;3,2;4,2;4,2];
[aBH,bBH]=find(blocks_list==1);
[aNH,bNH]=find(blocks_list==2);
[aBL,bBL]=find(blocks_list==3);
[aNL,bNL]=find(blocks_list==4);
for i=1:length(aBH)
    blocks_events(aBH(i),bBH(i),:)=BH_events(i,:);
    blocks_events(aNH(i),bNH(i),:)=NH_events(i,:);
    blocks_events(aBL(i),bBL(i),:)=BL_events(i,:);
    blocks_events(aNL(i),bNL(i),:)=NL_events(i,:);
end
%%

nsimsubs=500;
w=rand(nsimsubs,1)*3.3-2.3;
w=exp(w);
leak=rand(nsimsubs,1)*0.5;
noise=rand(nsimsubs,1)*3.3-2.3;
noise=exp(noise);
lapse=rand(nsimsubs,1)*0.2;

paraset=[w,leak,noise,lapse];
%figure;plot(opt1_events);hold on;plot(opt2_events);legend('opt1','opt2')
%%
nabandon=5;
for i=1:size(paraset,1)
 prob1(i,:,:)=selective_gating(squeeze(blocks_events(:,1,:)),squeeze(blocks_events(:,2,:)),paraset(i,1),paraset(i,2),paraset(i,3),paraset(i,4),nabandon);
end

%% generate choices
for i=1:size(prob1,1)
    for j=1:size(prob1,2)
        for k=1:size(prob1,3)
          choices(i,j,k)=randsample([1 0],1,1,[prob1(i,j,k),1- prob1(i,j,k)]);       
        end
    end
end
%%
save('simulated_data_SI_model.mat','choices','paraset','blocks_events')
