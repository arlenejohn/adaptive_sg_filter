function [GUE_MSE]=Find_risk_win(H,A,winlength,sigma,noisy_sig,numiter,len)
% winlength+1 is the total window length
% sigma is the standard deviation estimate
% noisy_sig is the noisy observation
% len is the length of the signal

th=-(winlength-1)/2:(winlength-1)/2;
GUE_MSE=zeros(1,len);

for m=1:len
    d=min(numiter,len-1)+m+th(1);
    e=min(numiter,len-1)+m+th(end);
    D=A*H;
    z=D*noisy_sig(d:e)';
    term1=z'*z;
    term2=-2*z'*noisy_sig(d:e)';
    term3=2*trace(D')*sigma^(2);
    term4=sum(noisy_sig(d:e).^2);
    GUE_MSE(1,m)=((term1+term2+term3+term4)/(winlength+1))-sigma^(2);
          
end