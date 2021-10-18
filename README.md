# adaptive_sg_filter
Code for the article "Adaptive Savitzky-Golay Filtering in Non-Gaussian Noise"

Denoising_main.m contains a commented description on how to call the functions that execute variations of the order adaptive and window adaptive SG filters discussed in the paper.

G-FL: den_win is the function that executes the G-FL algorithm or the window length adaptive SG filter with fixed order. It calculates the risk estimate using the Find_risk_win function and the standard deviation estimate using the sigma_estimate function.

G-FL-R: den_win_reg is the function that executes the G-FL-R algorithm or the window length adaptive SG filter with fixed order and regularization. It calculates the risk estimate using the Find_risk_win_reg function and the standard deviation estimate using the sigma_estimate function.

G-O: den_ord is the function that executes the G-O algorithm or the order adaptive SG filter with fixed window length. It calculates the risk estimate using the Find_risk_ord function and the standard deviation estimate using the sigma_estimate function.

G-O-R: den_ord_reg is the function that executes the G-O-R algorithm or the order adaptive SG filter with fixed window length and regularization. It calculates the risk estimate using the Find_risk_ord_reg function and the standard deviation estimate using the sigma_estimate function.

The aami3am.mat file is a clean ecg signal and the aami3am_noisy_L.mat file is the noisy ecg signal with noisy laplacian noise added to it, which is included for demonstrative purposes.
