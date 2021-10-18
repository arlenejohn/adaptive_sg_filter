function [denoised,order,GUE_MSE] = den_ord(M,pmax,noisy_sig,type,sigma)
%% Denoising Using Optimum order SG filter

%%%% Input
     % M: half-window length
     % pmax: Maximum Order
     % noisy_sig: Noisy signal (1XN)
     % type: Noise type ('G': Gaussian, 'L': Laplacian, 'U': Uniform)
     % sigma: noise standard deviation
%%%% Output
    % denoised: Denoised output 
    % Order: Optimum order at each time instants 
    % GUE_MSE:  GUE-MSE corresponding to the order values 1 to pmax, at each time instants
    
%% Function dependencies %%%%% IMPORTANT %%%%% 
    %%% Find_risk_ord.m
    %%% sigma_estimate.m
    
%%    

if(isrow(noisy_sig)~=1)
    noisy_sig=noisy_sig';
end

pmin=1;
[~,h]=size(noisy_sig);
H=cell(1,pmax);
A=cell(1,pmax);
winlen=2*M+1;

GUE_MSE=zeros(pmax,h);
denoised=zeros(1,h);
order=zeros(1,h);

%%
for p=1:pmax %calculating A and H matrices for differnt orders
    for x=1:winlen
       A{p}(x,:)=(-(winlen-1)/2+(x-1)).^(0:p);
    end
    H{p}=(A{p}'*A{p})\A{p}';
end
coefficients= cell(1,pmax); %the coefficient matrix

samp_noisy=[noisy_sig(h:-1:2) noisy_sig noisy_sig(h-1:-1:1)];
%% sigma estimation
if nargin<5
    sigma = sigma_estimate(median(abs(noisy_sig(2:end)-noisy_sig(1:end-1))),type);
end
%% Coefficient (a) determination
i=1;
while(i<=h)
    while(pmin<=pmax)
       r=samp_noisy(h-1+i-M:h-1+i+M);
       a=H{pmin}*r';
       coefficients{pmin}(i,:)=a';
       pmin=pmin+1;
    end
    i=i+1;
   pmin=1;
end
%% risk estimation at each point of reconstruction 
for t=1:pmax
 GUE_MSE(t,:)=Find_risk_ord(H{t},A{t},winlen,sigma,samp_noisy,h);
end
%% choosing optimal order at each reconstruction instant and evaluating the signal estimate
for k=1:h
    [~,index]=min(GUE_MSE(:,k));
    denoised(k)=coefficients{index}(k,1);
    order(k)=index;
end



