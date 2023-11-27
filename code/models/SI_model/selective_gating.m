function [estY, estX, cumX] = selective_gating(Xa,Xb,w,leak,noise,lapse,nabandon,start)
%function [estY, estX, cumX] = Model_SelectiveIntegration(modparam,X)
% parameters = [w leak noise lapse]
% Based on Glickman 2018
% w: selctive geting
% leak:  integration leak
%%
if nargin<7 nabandon=0;end
if nargin<8 start=0.5; end
%%
ntrials = size(Xa,1);
nsamp   = size(Xa,2);
% %
% disp(['ntrials : ',num2str(ntrials)])
% disp(['nsamples in a trial : ',num2str(nsamp)])
% Selective weighting v1/v2
diffX   = Xb - Xa; % R - L
wXb     = sigmf(diffX,[w,0]);
wXa     = 1 - wXb;
Ia      = wXa .* Xa;
Ib    	= wXb .* Xb;

% Accumulation v2 (same results as v1)
Ya(:,1)=ones(1,size(Ia,1)).*start;
Yb(:,1)=ones(1,size(Ib,1)).*start;

for i=2:nsamp
    Ya(:,i)= (1-leak)*Ya(:,i-1)+Ia(:,i-1);
    Yb(:,i)= (1-leak)*Yb(:,i-1)+Ib(:,i-1);
end
% Decision

diffY   = Ya(:,nabandon+1:end) - Yb(:,nabandon+1:end);
estY    = lapse + (1 - lapse*2).*sigmf(diffY,[noise,0]);
estX    = cat(3,Ia,Ib);
cumX    = cat(3,Ya,Yb);

end