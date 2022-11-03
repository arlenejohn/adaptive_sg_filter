function clean_signal = robust(noisy_sig)

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
    
M=9;
p=2;

%%

if(isrow(noisy_sig)~=1)
    noisy_sig=noisy_sig';
end


[~,h]=size(noisy_sig);

noisy_samp=[noisy_sig(M+1:-1:2) noisy_sig noisy_sig(h-1:-1:h-(M))];

[A,H]=find_H(M,p);

samp=noisy_samp';
winlength=2*M+1;
denoised=zeros(1,h);
th=-(winlength-1)/2:(winlength-1)/2;
for m=1:h
    d=min(M+1,h-1)+m+th(1)-1;
    e=min(M+1,h-1)+m+th(end)-1;
    a=H*samp(d:e);
    check=1;
    p_prev=samp(d:e);
    p=A*a;
    while check==1
        if ((norm((p_prev-p),2)/norm(p_prev,2))<0.0001)
            break
        end
        diff=1./max(0.0000002,abs(samp(d:e)-p));
        D=diag(diff);
        a=(A'*D*A)\(A'*D*samp(d:e));
        p_prev=p;
        p=A*a;
    end
    denoised(m)=a(1);
end
denoised(isnan(denoised))=0;
clean_signal= denoised;



