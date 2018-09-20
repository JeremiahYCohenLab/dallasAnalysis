function optoIDintan(filename, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('SamplingFreq', 20000)
p.addParameter('Pulses', 10)
p.addParameter('Trains', 10)
p.addParameter('PulseWidth', 10)
p.addParameter('PulseFreq', 10)
p.addParameter('ResponseWindow', 30)
p.parse(varargin{:});

%get session info
[root, sep] = currComputer();
[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];

%specify and make directories
sortedDataPath = [root animalName sep sessionFolder sep 'sorted' sep p.Results.Session sep];
unsortedDataPath = [root animalName sep sessionFolder sep 'ephys' sep p.Results.Session sep];
saveDir = [root animalName sep sessionFolder sep 'figures'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%get sorted and raw ephys data
optoFiles = dir(fullfile(unsortedDataPath,'*.rhd'));
sortedFiles = dir(fullfile(sortedDataPath,'*.txt'));

%get list of channels with sorted units
for i = 1:length(sortedFiles)
    tmpInds = strfind(sortedFiles(i).name, '_');
    chanList(i) = str2double(sortedFiles(i).name(tmpInds(1)+1:tmpInds(2)-1));
end
chanListU = unique(chanList);

%combine traces and DI data from raw files, only from relvant channels
laser = [];
traces = [];
for i = 1:length(optoFiles)
    [digInTemp, tracesTemp, ~] = readIntan(unsortedDataPath, optoFiles(i).name);
    laser = [laser digInTemp(8,:)];
    traces = [traces tracesTemp(chanListU,:)];
end
laserOnInds = find(laser(1:end-1) == 0 &  laser(2:end) > 0) + 1;
laserOffInds = find(laser(1:end-1) > 0 &  laser(2:end) == 0) + 1;

%create timestamps
tSamp = 1/p.Results.SamplingFreq * 1e3; % time per sample in ms
ts_interp = 0:tSamp:tSamp*(length(traces));
tSamp_us = 1/p.Results.SamplingFreq * 1e6;
ts_interp_us = 0:tSamp_us:tSamp_us*(length(traces));

laserOn = ts_interp(laserOnInds); %laser time on in ms
laserOff = ts_interp(laserOffInds);

%set window and stim paramaters
tB = 500;
tA = 500;
pulseInds = (1:p.Results.Pulses:p.Results.Pulses*p.Results.Trains);
respWin = p.Results.ResponseWindow;
rasterLength = length(-1*tB:(p.Results.Pulses*(1000/p.Results.PulseFreq)+tA));
    

%% 

for i = 1:length(sortedFiles)
    [cellName, ~] = strtok(sortedFiles(i).name, '.');
    spikeTimesTmp = [load(strcat(sortedDataPath, sortedFiles(i).name))]';
    
    [spikeInds, ~] = ismember(ts_interp_us, spikeTimesTmp);     %convert to us to get correct/whole indices
    spikeInds = find(spikeInds == 1);
    for j = 1:length(spikeInds)
        traceSlopeTmp = fliplr(diff(traces(i,spikeInds(j)-40:spikeInds(j)))); %get slopes of trace in 2ms before spike peak, flip
        tmpInd = find((traceSlopeTmp(1:end-1) < 0 & traceSlopeTmp(2:end) > 0) | (traceSlopeTmp(1:end-1) > 0 & traceSlopeTmp(2:end) < 0), 1) + 1;
        spikeTimes(j) = spikeTimesTmp(j) - tmpInd*tSamp_us;
    end
    spikeTimes = spikeTimes/1000;           %convert to ms from us
    
    spikeRast = [];
    for j = 1:p.Results.Trains
        spikeRast{j} = spikeTimes(spikeTimes > (laserOn(pulseInds(j)) - tB) & spikeTimes < (laserOff(pulseInds(j)+9) + tA));
        spontSpikeRast{j} = spikeTimes((spikeTimes > (laserOn(pulseInds(j)) - tB) & spikeTimes < laserOn(pulseInds(j))) |...
            (spikeTimes > (laserOff(pulseInds(j)+9) + respWin) & spikeTimes < (laserOff(pulseInds(j)+9) + tA)));
        if ~isempty(spikeRast{j})
            spikeRast{j} = spikeRast{j} - laserOn(pulseInds(j)); %puts in time relative to first light pulse
        end
    end
    
    %find times when there is no light for control comparison
    laserSham = [linspace(-tB, 0-1000/p.Results.PulseFreq, p.Results.Pulses/2) ...
        linspace(p.Results.Pulses*1000/p.Results.PulseFreq, rasterLength-1000/p.Results.PulseFreq, p.Results.Pulses/2)];
    
    spikeLat = [];
    spikeLatSham = [];
    lightSpikeTimes = [];
    for j = 1:p.Results.Trains              %for all pulses in all trains, find spikes within the response window
        for k = 1:p.Results.Pulses
            spikeRespTemp = spikeTimes(spikeTimes > laserOn(pulseInds(j)+k-1) & ...
                spikeTimes < laserOn(pulseInds(j)+k-1) + respWin);
            spikeRespTempSham = spikeTimes(spikeTimes > laserSham(k) + laserOn(pulseInds(j)) & ...
                spikeTimes < laserSham(k) + laserOn(pulseInds(j)) + respWin);
            if ~isempty(spikeRespTemp)
                spikeLat(j,k) = spikeRespTemp(1) - laserOn(pulseInds(j)+k-1);
                lightSpikeTimes = [lightSpikeTimes spikeRespTemp(1)];
            else
                spikeLat(j,k) = NaN;
            end
             if ~isempty(spikeRespTempSham)
                spikeLatSham(j,k) = spikeRespTempSham(1) - (laserSham(k) + laserOn(pulseInds(j)));
            else
                spikeLatSham(j,k) = NaN;
            end           
        end
    end
    avgSpikeLat = nanmean(spikeLat);    avgSpikeLatSham = nanmean(spikeLatSham);        %find average spikeLat and P(spike)
    spikeProb = mean(~isnan(spikeLat)); spikeProbSham = mean(~isnan(spikeLatSham));
    
    %using spike times, get indices for extracting spike waveforms
    [lightSpikeInds, ~] = ismember(ts_interp_us, lightSpikeTimes*1000);     %convert to us to get correct/whole indices
    lightSpikeInds = find(lightSpikeInds == 1);
    for j = 1:length(lightSpikeInds)
        lightSpikeMat(j,:) = [(lightSpikeInds(j)-10) : (lightSpikeInds(j)+20)];     %create mat for waveforms
    end
    
    spontSpikeTimes = [];
    for j = 1:length(spontSpikeRast)
        if ~isempty(spontSpikeRast{j})
            spontSpikeTimes = [spontSpikeTimes spontSpikeRast{j}];
        end
    end
    [spontSpikeInds, ~] = ismember(ts_interp_us, spontSpikeTimes*1000);
    spontSpikeInds = find(spontSpikeInds == 1);
    
    for j = 1:length(spontSpikeInds)
        spontSpikeMat(j,:) = [(spontSpikeInds(j)-10) : (spontSpikeInds(j)+20)];
    end
    
    cellTrace = traces(i,:);
    lightSpikeTraces = cellTrace(lightSpikeMat);       %extract waveforms from trace
    spontFlag = exist('spontSpikeMat', 'var');
    if spontFlag
        spontSpikeTraces = cellTrace(spontSpikeMat);
    end
    
    
    %% plot everything
    
    rasters = figure; subplot(4,3,[1:6]); hold on; title(strcat(filename, '_', cellName),'Interpreter','none')
    xlabel('Time (ms)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    x = linspace(0, ((p.Results.Pulses-1)*1000/p.Results.PulseFreq), p.Results.Pulses);
    xx = x + p.Results.PulseWidth;
    for j = 1:length(x)
        plotShaded([x(j) xx(j)],[0 0; 1+p.Results.Trains 1+p.Results.Trains],'b');
    end
    
    subplot(4,3,[7]); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 30]); xlim([0 p.Results.Pulses+1])
    plot(avgSpikeLat, 'b');
    plot(avgSpikeLatSham, 'k');
    legend('laser','control');
    
    subplot(4,3,10); hold on;
    xlabel('Pulse'); ylabel('P(spike)'); ylim([-0.1 1.1]); xlim([0 p.Results.Pulses+1])
    plot(spikeProb, 'b');
    plot(spikeProbSham, 'k');
    
    subplot(4,3,8); hold on;
    ylabel('Amplitude (uV)');
    plot([ts_interp(1:size(lightSpikeTraces, 2)) - .5], mean(lightSpikeTraces), '-b', 'linewidth', 2)
    
    if spontFlag
        subplot(4,3,11); hold on;
        xlabel('Times (ms)'); ylabel('Amplitude (uV)');
        plot([ts_interp(1:size(spontSpikeTraces, 2)) - .5], mean(spontSpikeTraces), '-k', 'linewidth', 2)
    end
    
    subplot(4,3,[9 12]); hold on
    xlabel('Time (ms)'); ylabel('Amplitude (uV)')
    errorfill([ts_interp(1:size(lightSpikeTraces, 2)) - .5], mean(lightSpikeTraces), std(lightSpikeTraces), 'b');
    if spontFlag
        errorfill([ts_interp(1:size(spontSpikeTraces, 2)) - .5], mean(spontSpikeTraces), std(spontSpikeTraces), 'k');
    end
 
    set(rasters, 'Position', get(0,'Screensize'))
    saveFigurePDF(rasters,[saveDir sep filename '_' cellName '_' p.Results.Session 'ID'])
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    