function [GUE_MSE]=Find_risk_ord_reg(H,A,winlength,sigma,noisy_sig,len,lambda)
% winlength - the total window length
% sigma - the standard deviation estimate
% samp - Observation
% len - is the length of the signal

if nargin<7
  lambda=12;  
end

th=-(winlength-1)/2:(winlength-1)/2;

GUE_MSE=zeros(1,len);

for m=1:len
    d=len-1+m+th(1);
    e=len-1+m+th(end);
    
    D=A*H;
    z=D*noisy_sig(d:e)';
    
    %term1=noisy_sig(d:e)*H'*A'*A*H*noisy_sig(d:e)';
    term1=z'*z;
 
    %term2=-2*noisy_sig(d:e)*H'*A'*((noisy_sig(d:e)'));
    term2=-2*z'*noisy_sig(d:e)';
    
    %term3=2*trace(H'*A')*sigma^(2);
    term3=2*trace(D')*sigma^(2);
    
    %term4=lambda*sum(diag(H'*A').^(2))*sigma^2;
    term4=lambda*sum(diag(D').^(2))*sigma^2;
    
    term5=sum(noisy_sig(d:e).^2);
    
    GUE_MSE(1,m)=((term1+term2+term3+term4+term5)/(winlength+1))-sigma^2;
        
end


