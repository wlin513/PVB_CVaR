%% checking behavior data quality
%online
%b1_bhbl;
%b2_nhnl;
%b3_bhnl;
%b4_nhbl;
%b5_bhnh;
%b6_blnl;
%b7_bihnh;
%b8_bilnl;
for i=1:size(data,2)
    %probabilities of choosing the same shape for each block
    tmp_optionchocies=transpose(reshape(data(i).choice,30,length(data(i).blockN)));
    k=1;
    for j=data(i).blockN
    opt1chocies(i,j,:)=tmp_optionchocies(k,:);
    k=k+1;
    end
    %probabilities of choosing the same side for each block
    tmp_sidechocies=mean(reshape(data(i).chosenside,30,length(data(i).blockN),1));
    k=1;
    for j=data(i).blockN
    sidechocies(i,j)=tmp_sidechocies(k);
    k=k+1;
    end
    %accuracies for each block
    tmp_optionchocies=transpose(reshape(data(i).choice,30,length(data(i).blockN)));
    tmp_correctopt=2-transpose(reshape(data(i).correctopt,30,length(data(i).blockN)));
    tmp_accuracies=tmp_correctopt==tmp_optionchocies;
    tmp_accuracies=double(tmp_accuracies);
    tmp_accuracies(tmp_correctopt==-1)=0.5;
    k=1;
    for j=data(i).blockN
    accuracy(i,j,:)=tmp_accuracies(k,:);
    k=k+1;
    end
    %rt for each block
    tmp_rt=transpose(reshape(data(i).rt,30,length(data(i).blockN)));
    k=1;
    for j=data(i).blockN
    rt(i,j,:)=tmp_rt(k,:);
    k=k+1;
    end
end
%calculating win stay and loss switch probilities and over-all switch
%probabilities
tmp1=accuracy(:,:,1:29);
tmpc2=opt1chocies(:,:,2:30);
tmpc1=opt1chocies(:,:,1:29);
tmp_ws=tmp1==1&tmpc2==tmpc1;
tmp_ls=tmp1==0&tmpc2~=tmpc1;
winstayprob=sum(tmp_ws,3)./sum(accuracy,3);
lossswitchprob=sum(tmp_ls,3)./sum(~accuracy,3);
wsls=winstayprob-lossswitchprob;
switchprobs=tmpc2~=tmpc1;
switchprobabilities=mean(switchprobs,2);
%
optionchocies=squeeze(mean(opt1chocies,3));
figure;imagesc(optionchocies>0.95|optionchocies<0.05)
imagesc(optionchocies==0)
imagesc(optionchocies==1)
tmp=(optionchocies==1|optionchocies==0);%(optionchocies>0.95|optionchocies<0.05);
find(sum(tmp(:,5:8),2)>0)
option_include=sum(tmp(:,5:8),2)==0;


% find(sum(tmp(:,5:8),2)>1)
% option_include=sum(tmp(:,5:8),2)<2;

imagesc(sidechocies==1|sidechocies==0)
find(sum(sidechocies==1|sidechocies==0,2)>0)
side_include=sum(sidechocies==1|sidechocies==0,2)==0;


accuracies=mean(accuracy,3);
figure;imagesc(accuracies)
accuracy_include=mean(accuracies(:,1:4),2)>0.6;

% for i=1:length(aa)
% totm(i)=data(aa(i)).total_money(end);
% end

    figure;
    for i=1:8
        scatter(ones(size(accuracies,1),1)*i,accuracies(:,i),[],'k','filled')
        hold on
    end
    
    plot(accuracies')
 %plot(mean(accuracies(:,1:8)))
 
% accuracy_include=accuracies(:,5)>=0.5;
% imagesc(accuracies(:,5)<0.5)
% 
% figure;
% bar(mean(accuracies,1))
% err=std(accuracies)/sqrt(size(subinfo,1));
% hold on
% er=errorbar(mean(accuracies,1),err);
% er.LineStyle='none';

f=figure;
bar(mean(optionchocies,1));
err=std(optionchocies)/sqrt(size(data,2));
hold on
er=errorbar(mean(optionchocies,1),err);
er.LineStyle='none';
ylabel('precentages of choosing option A');
%xticks([1:8])
xticklabels({'BHBL','NHNL','BHNL','NHBL','BHNH','BLNL','BiHNH','BiLNL'});
%saveas(f,'perc_ch_optA_beijing_online.png')
%% check RT differences between blocks
f=figure;
bar(mean(median(rt,3),1));
err=std(median(rt,3))/sqrt(size(data,2));
hold on
er=errorbar(mean(median(rt,3),1),err);
er.LineStyle='none';
ylabel('RT');
%xticks([1:8])
xticklabels({'BHBL','NHNL','BHNL','NHBL','BHNH','BLNL','BiHNH','BiLNL'});

aa=median(rt,3);
[h,p]=ttest(mean(aa(:,1:4),2),mean(aa(:,5:8),2))