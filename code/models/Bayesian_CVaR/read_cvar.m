
function [cvar1,cvar2]=read_cvar(eta,opt1_expv,opt2_expv)
    if eta<=0
        prc=(eta+1)*100;
        var1=prctile(opt1_expv,prc);
        for i=1:size(var1,2)
            for j=1:size(var1,3)
                cvar1(i,j)=mean(opt1_expv(opt1_expv(:,i,j)<=var1(1,i,j),i,j));
            end
        end
        
        var2=prctile(opt2_expv,prc);
        for i=1:size(var2,2)
            for j=1:size(var2,3)
                cvar2(i,j)=mean(opt2_expv(opt2_expv(:,i,j)<=var2(1,i,j),i,j));
            end
        end
    else
        prc=eta*100;
        var1=prctile(opt1_expv,prc);
        for i=1:size(var1,2)
            for j=1:size(var1,3)
                cvar1(i,j)=mean(opt1_expv(opt1_expv(:,i,j)>=var1(1,i,j),i,j));
            end
        end
        
        var2=prctile(opt2_expv,prc);
        for i=1:size(var2,2)
            for j=1:size(var2,3)
                cvar2(i,j)=mean(opt2_expv(opt2_expv(:,i,j)>=var2(1,i,j),i,j));
            end
        end
    end
end
