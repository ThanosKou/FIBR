function [throughput] = distBasedThroughput(OnPeriods, BS_distances)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

choose = @(samples) samples(randi(numel(samples)));
random_BS_dist = zeros(1,length(OnPeriods));

for i=1:length(OnPeriods)
    random_BS_dist(i) = choose(BS_distances);
end 


d0 = 1;

fc = 28;
n = 4.6;
sigma = 8.54;
FSPL = 32.4 + 20*log10(fc);

beta = 0.85; % part of bandwidth for transmission
gamma = 0.8;
W = 100*10^6;
k_B = 1.381*10^(-23);
T = 298; % temperature
noise_JN = 10*log10(W*k_B*T); %johnson nyquist

X_sigma = sigma.*randn(1,length(OnPeriods));

PL = FSPL + 10*n*log10(random_BS_dist/d0) + X_sigma;

signal = 27 ; % initial transmitted in dBm
signal = signal - PL; %dBm
signal = signal + 2*27; %beamforming gain, still in dB


SNR_dB = signal - noise_JN - 30; % to make it in dB
SNR = 10.^(SNR_dB/10);


capacity  = beta*W*log2(1+SNR);

simTime = 4*60*60;
throughput = sum(capacity.*OnPeriods)/simTime;

end 
