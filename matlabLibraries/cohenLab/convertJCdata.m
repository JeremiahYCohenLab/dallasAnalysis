% convertJCdata

cd('C:\Users\naoshigeuchida\Documents\MATLAB\Cohen')
load('brAll.mat')

NN = 95;    % number of neurons
TN = 234;   % max number of trials
SN = 559;   % max number of spikes

for i = 1:NN
    ts = [];
    OdorOn = [];
    ntrials = 0;
    for k = 1:TN
        for l = 1:SN
            if ~isnan(brAll(i,k,l))
                TS = 10*1000*k;
                ts = [ts (brAll(i,k,l)+TS)*10];   % 10^-4 sec
                OdorOn(k) = TS/1000;            % sec
                ntrials = k;
            end
        end
    end
    
    if k<10; str = ['00' num2str(i)];
    elseif k<100; str = ['0' num2str(i)];
    else str = num2str(i);
    end
    
    fn = ['T' str '.mat'];
    save(fn,'ts','OdorOn','ntrials')
end

