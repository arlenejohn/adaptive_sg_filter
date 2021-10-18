function[x,sigma]=add_noise_2(s,SNR_dB,noise)
% adds noise of desired SNR to the signal of desired ditribution
 %% input
 %   s-clean signal
 %  SNR_dN - input SNR in dB
 %  noise - type of noise distribution ('G','L','U')
 %% output
 % x - noisy signal
 %% 
squares=(norm(s).^(2))/length(s);
Variance=squares/(10^(SNR_dB/10));
if(noise=='G')
    nois=randn(1,length(s));
    x=s+sqrt(Variance)*nois;
    sigma=sqrt(Variance);
end
if(noise=='L')
    nois=laprnd(1,length(s), 0, 1);
    x=s+sqrt(Variance)*nois;
    sigma=sqrt(Variance);
end
if(noise=='U')
    a=-0.5;
    b=0.5;
    nois = a + (b-a).*rand(1,length(s));
    x=s+sqrt(Variance)*(sqrt(12)*nois);
    sigma=sqrt(Variance);
end