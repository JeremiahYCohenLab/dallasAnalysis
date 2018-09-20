% This function computes a spike density function
% x is a spike train of zeroes and ones
% Filter is based on:
%   Thompson KG, Hanes DP, Bichot NP, Schall JD. Perceptual and motor
%     processing stages identified in the activity of macaque frontal eye
%     field neurons during visual search. J Neurophysiol 76:
%     4040-4055, 1996.
% Written by Jeremiah Y. Cohen

function mySDF = sdf(x)
    % This is a postsynaptic-potential-like filter to convolve with spikes
    myPSP = (1 - exp(-(1 : 80))) .* (exp(-(1 : 80)/20));
    convf = 1000 * myPSP./sum(myPSP);
    
    % Convolve the filter with the spike train
    mySDF = conv(convf, x);
