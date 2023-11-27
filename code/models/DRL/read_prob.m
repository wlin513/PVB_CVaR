function prob = read_prob(dist1,dist2,nsample)
%read in two distributions and calculate the probability that dist1 higher
%than dist2 using permutation
if nargin<3
    nsample=length(dist1)^2;
end

 if std(dist1)==0 && std(dist2)==0 && mean(dist1)==mean(dist2)
     prob=0.5;
 else  
    sample1=randi([1,length(dist1)],nsample,1);
    sample2=randi([1,length(dist2)],nsample,1);
    dist1_samples = dist1(sample1);
    dist2_samples = dist2(sample2);
    prob=sum(dist1_samples>dist2_samples)/nsample;
 end
end

