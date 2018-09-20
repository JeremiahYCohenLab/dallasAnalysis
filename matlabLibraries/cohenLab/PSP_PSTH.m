function [mySDF, t] = PSP_PSTH(TS,trigger,trials,window)

TS = TS * 10^-4;

plt = 'n';
err = 2;
N = length(trials);

for i = 1:length(trials)
    tt = trials(i);
    t0 = window(1)-1;   % take data from + or - 1 sec
    t1 = window(2)+1;
    sp = TS-(trigger(tt));
    data(i).times = (sp(sp>=t0 & sp<=t1))';
end

n = 1000*(window(2)-window(1));  % ms
x = zeros(1:length(trials),n);
for i = 1:N
    d = data(i).times-(window(1));
    d = d*1000;
    d = d(d<(window(2)-window(1))*1000 & d>0);
    d = round(d);
    x(i,d) = 1;
end

t = window(1):0.001:window(2);  % 1 ms-step


% % PSP filter
% n = 80;
% tau = 20;
% f = (1 - exp(-(1 : n))) .* (exp(-(1 : n)/tau));
% ind = 1:length(t);

% Gaussian filter
sigma = 20;
xf = -4*sigma:4*sigma;
f = normpdf(xf,0,sigma);
n = length(f);
ind = [1:length(t)] + round(n/2);

mySDF = [];
convf = 1000 * f./sum(f);
for i = 1:N
% Convolve the filter with the spike train
    mySDF(i,:) = conv(convf, x(i,:));
end

mySDF = mySDF(:,ind);

