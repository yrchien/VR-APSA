function y = BG_noise(p,sigma,GINR,N)

% p : probability of impulsive noise occurrence 
% sigma : standard deviation of AWGN
% GINR : Gaussian to Impulsive noise ratio
% N : length of noise

b = binornd(1,p,1,N);
g = randn(1,N)*sigma*(1/GINR);
y = b.*g;