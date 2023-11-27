function [v1,v2,s1,s2]=PEIRS(opt1_events,opt2_events,alphaQ,alphaS,startQ,startS)

%

v1=zeros(size(opt1_events));
pe1=zeros(size(opt1_events));
s1=zeros(size(opt1_events));
v1(:,1)=startQ;
s1(:,1)=startS;
v2=zeros(size(opt2_events));
pe2=zeros(size(opt2_events));
s2=zeros(size(opt1_events));
v2(:,1)=startQ;
s2(:,1)=startS;
for b=1:size(opt1_events,1)
    for i=1:size(opt1_events,2)-1
      pe1(b,i)=opt1_events(b,i)-v1(b,i);
      v1(b,i+1)=v1(b,i)+alphaQ*pe1(b,i);
      s1(b,i+1)=s1(b,i)+alphaS*(abs(pe1(b,i))-s1(b,i));
      
      pe2(b,i)=opt2_events(b,i)-v2(b,i);
      v2(b,i+1)=v2(b,i)+alphaQ*pe2(b,i);
      s2(b,i+1)=s2(b,i)+alphaS*(abs(pe2(b,i))-s2(b,i));     
    end    
end
