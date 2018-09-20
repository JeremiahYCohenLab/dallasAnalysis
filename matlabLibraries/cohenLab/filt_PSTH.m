function [mySDF, t] = filt_PSTH(TS,trigger,trials,window,ftype,param)

% filter type:
% Gaussian: param: sigma
% PSP: param: tau

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

% filter
if strncmp(ftype,'G',1)  % Gaussian filter
    sigma = param; % sigma
    xf = -4*sigma:4*sigma;
    f = normpdf(xf,0,sigma);
    n = length(f);
    ind = [1:length(t)] + round(n/2);
    
elseif strncmp(ftype,'B',1)  % Gaussian filter
    W = param; % window width (ms)
    f = ones(1,round(W));
    n = length(f);
    ind = [1:length(t)] + round(n/2);
    
elseif strncmp(ftype,'P',1) % PSP filter
    n = 80;
    tau = 20;
    f = (1 - exp(-(1 : n))) .* (exp(-(1 : n)/tau));
    ind = 1:length(t);
    
else
    disp('Choose Gaussian or PSP')
end


mySDF = [];
convf = 1000 * f./sum(f);
for i = 1:N
% Convolve the filter with the spike train
    mySDF(i,:) = conv(convf, x(i,:));
end

mySDF = mySDF(:,ind);

