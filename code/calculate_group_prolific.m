%% 
clear all;
clc;
%%
addpath('models\PNPE')
addpath('models\1LR1b')
addpath('models\SI_model')
addpath('models\DRL')
addpath('models\CVaR')
addpath('models\DRL_LR')
addpath('models\PEIRS\')
addpath('utility\')
addpath('models\Bayesian_CVaR\')
addpath('models\Bayesian_CVaR_eta_beta\')
addpath('models\Bayesian_CVaR_eta_priorvar\')
addpath('models\concave_UTIL\')
addpath('models\convex_UTIL\')
addpath('models\inverse_s_shape_UTIL\')
%cd('D:\2021_DistRL\code');
getfolders;
%% read in behavior and questionnaire data
dem_info=readtable([datadir,'prolific_export_prolific_beijing.csv']);

sublistfile=dir([datadir,'\prolific_beijing_version\*.csv']);
for i=1:size(sublistfile,1)
    i
    dfile=[sublistfile(i).folder,'\',sublistfile(i).name];
    [data(i),DAS(i),SDS(i),RRS(i),tSTAI(i),IUS(i)]=read_csv_prolific_beijing(dfile);
    tmp=readtable(dfile);
    pid=tmp.PROLIFIC_PID{1};
    j=find(strcmp(dem_info.ParticipantId,pid));
    demographic(i).age=dem_info.Age(j);
    demographic(i).ethnicity=dem_info.EthnicitySimplified(j);
    demographic(i).sex=dem_info.Sex(j);
end
%% calculate questionnaire scores and correlate them with model estimates

%mark where SDS should be reversed
sds_rev=zeros(20,1);
for i=[2,5,6,11,12,14,16,17,18,20]
    sds_rev(i,1)=1;
end
%mark where SDS should be reversed
tstai_rev=zeros(20,1);
for i=[1,3,6,7,10,13,14,16,19]
    tstai_rev(i,1)=1;
end
ques_include=ones(size(data,2),1);
for i=1:size(data,2)
    if sum(isnan(DAS(i).raw))==0
        data(i).DAS=sum(DAS(i).raw);
    end
    if sum(isnan(RRS(i).raw))==0
        data(i).RRS=sum(RRS(i).raw);    
    end
    if sum(isnan(SDS(i).raw))==0
        tmp=SDS(i).raw;
        if length(unique(tmp))==1
            ques_include(i)=0;
        end
        tmp(sds_rev==1)=5-tmp(sds_rev==1);
        data(i).SDS=sum(tmp);     
    end
    if sum(isnan(tSTAI(i).raw))==0
        tmp=tSTAI(i).raw;
        if length(unique(tmp))==1
            ques_include(i)=0;
        end
        tmp(tstai_rev==1)=5-tmp(tstai_rev==1);
        data(i).tSTAI=sum(tmp);         
    end
    
%DAS(i).total=sum(DAS(i).raw);
IUS_f1(i,1)=sum(IUS(i).raw([1, 2, 3, 9, 12, 13, 14, 15, 16, 17, 20, 22, 23, 24,  25],1));
IUS_f2(i,1)=sum(IUS(i).raw([4, 5, 6, 7, 8, 10, 11, 18, 19, 21, 26,  27],1));
data(i).IUS=IUS_f1(i,1)+IUS_f2(i,1);
end
ques=struct2table(data);
%%
check_beh_data_quality_prolific
save([datadir,'data_prolific_beijing_version.mat'],'data','side_include','option_include','accuracy_include','ques_include','optionchocies')
save('C:\Users\wlin\OneDrive - University College London\2021_DistRL\analyses_beijing_beh_data\data\experiment\data_prolific_beijing_version.mat','data','side_include','option_include','accuracy_include','ques_include','optionchocies')

%%
data_include=option_include&accuracy_include;


plotcorr(mean(optionchocies(data_include,5:8),2),ques.tSTAI(data_include),'pro-v','trait anxiety')
plotcorr(mean(optionchocies(data_include,5:8),2),ques.DAS(data_include),'pro-v','DAS')
plotcorr(mean(optionchocies(data_include,5:8),2),ques.RRS(data_include),'pro-v','RRS')
plotcorr(mean(optionchocies(data_include,5:8),2),ques.SDS(data_include),'pro-v','SDS')
plotcorr(mean(optionchocies(data_include,5:8),2),ques.IUS(data_include),'pro-v','IUS')
f3=figure;
quesallz=zscore(ques.RRS(data_include))+zscore(ques.tSTAI(data_include))+zscore(ques.SDS(data_include));
plotcorr(quesallz,mean(optionchocies(data_include,5:8),2),1,1,'sum zscores of all questionnaires','probabilities of choosing the broader options');


plotcorr(mean(optionchocies(:,5:8),2),ques.tSTAI,'pro-v','trait anxiety')
plotcorr(mean(optionchocies(:,5:8),2),ques.DAS,'pro-v','DAS')
plotcorr(mean(optionchocies(:,5:8),2),ques.RRS,'pro-v','RRS')
plotcorr(mean(optionchocies(:,5:8),2),ques.SDS,'pro-v','SDS')
plotcorr(mean(optionchocies(:,5:8),2),ques.IUS,'pro-v','IUS')

plotcorr(ques.tSTAI,ques.RRS,'trait anxiety','RRS')
plotcorr(ques.SDS,ques.RRS,'SDS','RRS')
%plot risk-seeking for both-high and both-low
figure;scatter(mean(optionchocies(data_include,[5,7]),2),mean(optionchocies(data_include,[6,8]),2));hold on;fplot(@(x) x,[0 1])

