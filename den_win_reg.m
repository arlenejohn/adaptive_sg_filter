function [denoised,window,GUE_MSE] = den_win_reg(Mmax,p,noisy_sig,type,lamda,sigma)

%% Denoising Using Optimum length SG filter with regularizer

%%%% Input
     % Mmax: Maximum half-window length
     % p:  Order
     % noisy_sig: Noisy signal (1XN)
     % type: Noise type ('G': Gaussian, 'L': Laplacian, 'U': Uniform)
     % lamda: Regularization parameter
     % sigma: noise standard deviation
%%%% Output
    % denoised: Denoised output 
    % window: Optimum window at each time instants 
    % GUE_MSE: Regularized GUE-MSE corresponding to Mmin to Mmax, at each time instants
    
%% Function dependencies %%%%% IMPORTANT %%%%% 
    %%% Find_risk_win_reg.m
    %%% sigma_estimate.m
    
Mmin=2; %% minimum windowlength 

%%

if(isrow(noisy_sig)~=1)
    noisy_sig=noisy_sig';
end

ow=Mmin:Mmax;
ow_srch_len=length(ow);

[~,h]=size(noisy_sig);
denoised=zeros(1,h);
window=zeros(1,h);

GUE_MSE=zeros(ow_srch_len,h);
H=cell(1,ow_srch_len);
A=cell(1,ow_srch_len);
coefficients=cell(1,ow_srch_len);

noisy_samp=[noisy_sig(Mmax+1:-1:2) noisy_sig noisy_sig(h-1:-1:h-(Mmax))];

%% sigma estimation
if nargin<6
    sigma = sigma_estimate(median(abs(noisy_sig(2:end)-noisy_sig(1:end-1))),type);
end

if nargin<5
    lamda=12;
end

%% calculating A and H matrices for different filter lengths
c=1;
for k=Mmin:Mmax
   [A{c}, H{c}]=find_H(k,p);
   c=c+1;
end
%% Coefficient (a) determination
i=1;
M_iter=Mmin;
c=1;
while(i<=h)
    while(M_iter<=Mmax)
       r=noisy_samp(min(Mmax,h-1)+i-M_iter:min(Mmax,h-1)+i+M_iter);
       a=H{c}*r';
       coefficients{c}(i,:)=a';
       M_iter=M_iter+1;
       c=c+1;
    end
    i=i+1;
    M_iter=Mmin;
    c=1;
end
%% risk estimation at each point of reconstruction 
j=1;
 for t=Mmin:Mmax
  GUE_MSE(j,:)=Find_risk_win_reg(H{j},A{j},2*t+1,sigma,noisy_samp,Mmax,h,lamda);
 j=j+1;
 end

 %% choosing optimal filter length at each reconstruction instant and evaluating the signal estimate
for k=1:h
    [~,index]=min(GUE_MSE(:,k));
    denoised(k)=coefficients{index}(k,1);
    window(k)=ow(index);
end


    