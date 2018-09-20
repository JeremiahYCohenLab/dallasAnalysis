function [spikes,ntrials,index] = raster12(TS,Trigger,trials,window,color,n,linewidth)

%RASTER12    - raster plot
%
% [spikes,ntrials,index] = raster12(TS,TrialStart,trials,window,color,n,linewidth)
%
%   input
%     - TS:         time stamps from spike data
%     - TrialStart: trial start (sec)
%     - trigger:    e.g. TrialStart, OdorPokeIn, OdorValveOn, OdorPokeOut, WaterPokeIn, WaterPokeOut etc.
%     - trials:     trial numbers
%     - window:     time wondow from the trigger (e.g. [-2 2])
%     - color
%     - N:          the starting y value for plotting
%     - linewidth
%     
%   output
%

if nargin < 7; linewidth = 0.5; end
if nargin < 6; n = 1; end
if nargin < 5; color = 'k'; end
    
h = 1;
spikes = [];
index = [];
k = n;
TS = TS * 10^(-4);
ntrials = 0;
for i = trials
    if ~isnan(Trigger(i))
        
        TRIG = Trigger(i);
        
        new_spikes = TS(find(TS>=(TRIG+window(1)) & TS<=(TRIG+window(2))))-TRIG;
        if ~isempty(new_spikes)
            if size(new_spikes,2) > 1; new_spikes = new_spikes'; end
            spikes = [spikes;new_spikes];
            index = [index;ones(length(new_spikes),1)*k];
        end
        k = k + 1;
        ntrials = ntrials + 1;
    end
end

h = line([spikes'; spikes'],[index'-h/2; index'+h/2],'color',color);
set(h,'linewidth',linewidth)