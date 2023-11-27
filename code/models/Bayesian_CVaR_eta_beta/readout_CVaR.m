function r=readout_CVaR(data,eta,step)
%data first dimension should be time, second should be distribution
%resolution (number of samples)
   if nargin<3
     step=1.01;
   end
    x=0.01:0.98/(size(data,2)-1):0.99;    
    
    for i=1:size(data,1)
        cdf=cumsum(data(i,:));
        %figure;plot(cdf)
        xp=nan;
        delta=0.001;
        if eta>0
            while isnan(xp)
                xp=round(median(find(cdf<(eta+delta)&cdf>(eta-delta))));
                delta=delta*step;
            end
            r(i)=sum(x(xp:end).*data(i,xp:end))/sum(data(i,xp:end));
        else
            if eta<0
                while isnan(xp)
                xp=round(median(find(cdf<(1+eta+delta)&cdf>(1+eta-delta))));
                delta=delta*step;
                end
                r(i)=sum(x(1:xp).*data(i,1:xp))/sum(data(i,1:xp));
            else
                r(i)=sum(data(i,:).*x);
            end  
        end
    end
end

