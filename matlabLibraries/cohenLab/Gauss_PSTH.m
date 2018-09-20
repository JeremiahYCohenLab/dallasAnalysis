function [R,t,E] = Gauss_PSTH(TS,trigger,trials,window,t_step,sigma)

TS = TS * 10^-4;

plt = 'n';
err = 2;
t = window(1):t_step:window(2);

for i = 1:length(trials)
    t = trials(i);
    t0 = window(1)-1;   % take data from + or - 1 sec
    t1 = window(2)+1;
    sp = TS-(trigger(t));
    data(i).times = (sp(sp>=t0 & sp<=t1))';
end

[R,t,E] = psth(data,sigma,plt,window,err); % Chronux toolbox