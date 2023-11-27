function rrdist=bayesian_update(data,priora,priorb,alpha,updatevar,N)
    %% Learning the values using bayesian rules - Fixed Belief Model(FBM)
    x=0.001:0.998/(N-1):0.999;
    x=x';
    prior=betapdf(x,priora,priorb);
    prior=prior/sum(prior);
    %figure;plot(prior)

    rdist=prior;
    rrdist(1,:)=prior;



    for i=2:length(data)


        rt=data(i-1);

        eventb=(rt - updatevar + rt .*updatevar - rt .^2*2 + rt .^3)./updatevar;
        eventa=-(rt .*(rt .^2 - rt  + updatevar))./updatevar;
        eventdist=betapdf(x,eventa,eventb)+1;

 
        %figure;plot(eventdist)
       
        rdist=alpha*rdist+(1-alpha)*prior;
        rdist=rdist.*eventdist;
        rdist=rdist/sum(rdist);
        %figure;plot(rdist)
        rrdist(i,:)=rdist;
    end
   
end