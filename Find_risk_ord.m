function  GUE_MSE=Find_risk_ord(H,A,winlength,sigma,noisy_sig,len)
%%%% Calculate unbiased estimate of MSE %%%%%

%%Input
    % winlength - the total window length
    % sigma - the standard deviation estimate
    % noisy_sig - the signal
    % len - the length of the signal

%% Output
       %%% GUE_MSE- Unbiased estimate of MSE
         
GUE_MSE=zeros(1,len);
th=-(winlength-1)/2:(winlength-1)/2;
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
        
    term4=sum(noisy_sig(d:e).^2);
    
    GUE_MSE(1,m)=((term1+term2+term3+term4)/(winlength+1))-sigma^(2);
           
end


