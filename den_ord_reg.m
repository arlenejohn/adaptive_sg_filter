function [denoised,order,GUE_MSE] = den_ord_reg(M,pmax,noisy_sig,type,lamda,sigma)
%% Denoising Using Optimum order SG filter

%%%% Input
     % M: half-window length
     % pmax: Maximum Order
     % noisy_sig: Noisy signal (1XN)
     % type: Noise type ('G': Gaussian, 'L': Laplacian, 'U': Uniform)
     % lamda: Regularization parameter
     % sigma: noise standard deviation
%%%% Output
    % denoised: Denoised output 
    % Order: Optimum order at each time instants 
    % GUE_MSE:  Regularized GUE-MSE corresponding to the order values pmin to pmax, at each time instants
    
%% Function dependencies %%%%% IMPORTANT %%%%% 
    %%% Find_risk_ord_reg.m
    %%% sigma_estimate.m

%%
pmin=1;

if(isrow(noisy_sig)~=1)
    noisy_sig=noisy_sig';
end

[~,h]=size(noisy_sig);
H=cell(1,pmax);
A=cell(1,pmax);
winlen=2*M+1;

GUE_MSE=zeros(pmax,h);
denoised=zeros(1,h);
order=zeros(1,h);

%% calculating A and H matrices for differnt orders
for p=1:pmax 
    for x=1:winlen
       A{p}(x,:)=(-(winlen-1)/2+(x-1)).^(0:p);
    end
    H{p}=(A{p}'*A{p})\A{p}';
end
coefficients= cell(1,pmax); %the coefficient matrix

noisy_samp=[noisy_sig(h:-1:2) noisy_sig noisy_sig(h-1:-1:1)];

%% sigma estimation
if nargin<6
    sigma = sigma_estimate(median(abs(noisy_sig(2:end)-noisy_sig(1:end-1))),type);
end
if nargin<5
    lamda=12;
end

%% Coefficient (a) determination
i=1;
while(i<=h)
    while(pmin<=pmax)
       r=noisy_samp(h-1+i-M:h-1+i+M);
       a=H{pmin}*r';
       coefficients{pmin}(i,:)=a';
       pmin=pmin+1;
    end
    i=i+1;
   pmin=1;
end
%% risk estimation at each point of reconstruction 
for t=1:pmax
 GUE_MSE(t,:)=Find_risk_ord_reg(H{t},A{t},winlen,sigma,noisy_samp,h,lamda);
end
%% choosing optimal order at each reconstruction instant and evaluating the signal estimate
for k=1:h
    [~,index]=min(GUE_MSE(:,k));
    denoised(k)=coefficients{index}(k,1);
    order(k)=index;
end



