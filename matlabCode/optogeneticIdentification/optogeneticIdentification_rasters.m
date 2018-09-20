function [latencyToResponse] = optogeneticIdentification_rasters(filename, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('filtFlag', 0)
p.addParameter('delayComp',0)
p.addParameter('TTdelay', 0)
p.addParameter('CSCdelay', 0)
p.addParameter('spikeDelay', 0)
p.addParameter('intanFlag', 0)
p.parse(varargin{:});


[root, sep] = currComputer();

[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'opto'];

if p.Results.filtFlag == 1
    filtDataPath = [root animalName sep sessionFolder sep 'ephys' sep 'opto' sep 'filtered' sep];
end

unsortedDataPath = [root animalName sep sessionFolder sep 'ephys' sep 'opto' sep];


if exist(strcat(sessionDataPath, sep, 'laserStruct.mat'),'file')
    load(strcat(sessionDataPath, sep, 'laserStruct.mat'))
else
    [laserStruct] = findLaserEvents(filename);
    load(strcat(sessionDataPath, sep, 'laserStruct.mat'))
end


if p.Results.delayComp == 1
    p.Results.spikeDelay = p.Results.TTdelay - p.Results.CSCdelay;
end

% find cluster data
spikeFields = fields(laserStruct);
clust = find(~cellfun(@isempty,strfind(spikeFields,'SS')) & ~cellfun(@isempty,strfind(spikeFields,'TT')));

sCon = 0.000000030518510385491027*1e6; %converts spike amplitude to uV
tB = 500; % ms before first pulse
tF = 500; % ms after last pulse
numPulses = 10; %number of pulses in each train
numTrains = length(laserStruct.laserOn)/numPulses;
rasterLength = length(-1*tB:laserStruct.laserOff(numPulses)-laserStruct.laserOn(1)+tF);

saveDir = [root animalName sep sessionFolder sep 'figures'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

for i = clust'
%     CSCchann = spikeFields{i}(strfind(spikeFields{i},'CSC'):end);
    
    TTchann = str2double(spikeFields{i}(3));
    CSCchann{1} = ['CSC' num2str(4*TTchann-3)];
    CSCchann{2} = ['CSC' num2str(4*TTchann-2)]; 
    CSCchann{3} = ['CSC' num2str(4*TTchann-1)];
    CSCchann{4} = ['CSC' num2str(4*TTchann-0)];
    trainInd = 1:numPulses:length(laserStruct.laserOn);
    
    latencyToResponse = []; % latency of light-evoked responses
    numResponses = []; % success of light-evoked response
    spikeTimes_light_TT = []; % light-evoked spike in TT time
    spikeTimes_light_CSC = []; % light-evoked spike in CSC time; always the first index of the 1x512 spike vector
    spikeCSC_light = []; % light-evoked spike waveform; 1x512 spike vector
    spikeInd_light = []; % index in the 1x512 spike vector that corresponds to spike

    spikeTimes_spont_TT = [];
    spikeTimes_spont_CSC = [];
    spikeCSC_spont = [];
    spikeInd_spont = [];
    
    for j = trainInd
        lL = laserStruct.laserOn(j) - tB;
        uL = laserStruct.laserOff(j+numPulses-1) + tF;
        spikeInds = laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i}) >= lL & laserStruct.(spikeFields{i}) < uL) - laserStruct.laserOn(j) + tB;
        spikeRast(ceil(j/numPulses),1) = {spikeInds};
        
        laserOnInds = laserStruct.laserOn(laserStruct.laserOn >= lL & laserStruct.laserOn < uL) - laserStruct.laserOn(j) + tB;
        laserOffInds = laserStruct.laserOff(laserStruct.laserOff >= lL & laserStruct.laserOff < uL) - laserStruct.laserOn(j) + tB;
        laserRast(ceil(j/numPulses),1) = {sort([laserOnInds laserOffInds])};
        
        for k = 1:numPulses
            tempSpike = laserStruct.(spikeFields{i})(find(laserStruct.(spikeFields{i}) > laserStruct.laserOn(j+k-1),1));
            tempLat = tempSpike - laserStruct.laserOn(j+k-1);
            if tempLat < 30 %30ms limit to evoke a spike
                latencyToResponse(ceil(j/numPulses),k) = tempLat - p.Results.TTdelay/1000;
                numResponses(ceil(j/numPulses),k) = 1;
%                 spikeTimes_light_TT(ceil(j/numPulses),k) = tempSpike; %(tempSpike*1000)-spikeDelay;
                spikeTimes_light_TT = [spikeTimes_light_TT tempSpike];
            else
                latencyToResponse(ceil(j/numPulses),k) = NaN;
                numResponses(ceil(j/numPulses),k) = 0;
%                 spikeTimes_light_TT(ceil(j/numPulses),k) = NaN;
%                 spikeTimes_light_TT = [spikeTimes_light_TT tempSpike];
            end
            
            if k < numPulses/2 + 1
                laserSham = lL + (k-1)*(tB/(numPulses/2));
            end
            if k > numPulses/2
                laserSham = uL - (k-1)*(tB/(numPulses/2));
            end
            tempSpikeSham = laserStruct.(spikeFields{i})(find(laserStruct.(spikeFields{i}) > laserSham, 1));
            tempLatSham = tempSpikeSham - laserSham;
            if tempLatSham < 30 %30ms limit to evoke a spike
                latencyToResponseSham(ceil(j/numPulses),k) = tempLatSham - p.Results.TTdelay/1000;
                numResponsesSham(ceil(j/numPulses),k) = 1;
            else
                latencyToResponseSham(ceil(j/numPulses),k) = NaN;
                numResponsesSham(ceil(j/numPulses),k) = 0;
            end
        end
        if j == 1
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i})  < lL)];
        elseif j == trainInd(end)
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i})  >= uL)];
        else
            nuL = laserStruct.laserOff(j+numPulses) - tB;
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i}) >= uL & laserStruct.(spikeFields{i}) < nuL)];
        end
    end
    spikeTimes_light_TT = 1000*spikeTimes_light_TT - p.Results.spikeDelay;
    
    if ~isempty(spikeTimes_light_TT)
        if p.Results.filtFlag == 1 
            if ~exist(strcat(filtDataPath,CSCchann{1},'.ncs'),'file')       %if filtered traces don't exit run function to filter CSC traces
                filterCSCtraces(unsortedDataPath, CSCchann, spikeTimes_light_TT);
                unsortedDataPath = filtDataPath;                            %change path to draw from filtered traces
            else
               unsortedDataPath = filtDataPath;
            end
        end
        
        [spikeTimes_light_CSC, spikeCSC_light] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{1},'.ncs'),[1 0 0 0 1],0,5,spikeTimes_light_TT);
        if length(spikeTimes_light_TT) > length(spikeTimes_light_CSC) % sometimes, spikeTimes_light_TT has indices that don't correspond to CSC data
            m = 1;
            while length(spikeTimes_light_TT) > length(spikeTimes_light_CSC)
                if spikeTimes_light_TT(m) - spikeTimes_light_CSC(m) < 0
                    spikeTimes_light_TT(m) = [];
                    m = m - 1;
                end
                m = m + 1;
            end
        end  
        
        [spikeCSC_light(:,:,2)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{2},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        [spikeCSC_light(:,:,3)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{3},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        [spikeCSC_light(:,:,4)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{4},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        spikeInd_light = floor((spikeTimes_light_TT - spikeTimes_light_CSC)/31.25);
       
    end

    spikeTimes_spont_TT = 1000*spikeTimes_spont_TT - p.Results.spikeDelay;
    [spikeTimes_spont_CSC, spikeCSC_spont] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{1},'.ncs'),[1 0 0 0 1],0,5,spikeTimes_spont_TT);
    if length(spikeTimes_spont_TT) > length(spikeTimes_spont_CSC) % sometimes, spikeTimes_spont_TT has indices that don't correspond to CSC data
        m = 1;
        while length(spikeTimes_spont_TT) > length(spikeTimes_spont_CSC)
            if spikeTimes_spont_TT(m) - spikeTimes_spont_CSC(m) < 0
                spikeTimes_spont_TT(m) = [];
                m = m - 1;
            end
            m = m + 1;
        end
    end
    [spikeCSC_spont(:,:,2)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{2},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_spont_TT);
    [spikeCSC_spont(:,:,3)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{3},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_spont_TT);
    [spikeCSC_spont(:,:,4)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{4},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_spont_TT);
    spikeInd_spont = floor((spikeTimes_spont_TT - spikeTimes_spont_CSC)/31.25);

    lightSpikeWaveforms = [];
    spontSpikeWaveforms = [];

    padding = 32; % indices to left and right of spike to record
    for l = 1:length(spikeInd_light)        
        if spikeInd_light(l) > padding && spikeInd_light(l) < 512-padding
            lightSpikeWaveforms(:,l,1) = spikeCSC_light(spikeInd_light(l) - padding:spikeInd_light(l) + padding, l, 1);
            lightSpikeWaveforms(:,l,2) = spikeCSC_light(spikeInd_light(l) - padding:spikeInd_light(l) + padding, l, 2);
            lightSpikeWaveforms(:,l,3) = spikeCSC_light(spikeInd_light(l) - padding:spikeInd_light(l) + padding, l, 3);
            lightSpikeWaveforms(:,l,4) = spikeCSC_light(spikeInd_light(l) - padding:spikeInd_light(l) + padding, l, 4);
        elseif spikeInd_light(l) <= padding
            lightSpikeWaveforms(:,l,1) = [NaN(padding-spikeInd_light(l)+1,1); spikeCSC_light(1:spikeInd_light(l)+padding, l, 1)];
            lightSpikeWaveforms(:,l,2) = [NaN(padding-spikeInd_light(l)+1,1); spikeCSC_light(1:spikeInd_light(l)+padding, l, 2)];
            lightSpikeWaveforms(:,l,3) = [NaN(padding-spikeInd_light(l)+1,1); spikeCSC_light(1:spikeInd_light(l)+padding, l, 3)];
            lightSpikeWaveforms(:,l,4) = [NaN(padding-spikeInd_light(l)+1,1); spikeCSC_light(1:spikeInd_light(l)+padding, l, 4)];
        elseif spikeInd_light(l) >= 512-padding
            lightSpikeWaveforms(:,l,1) = [spikeCSC_light(spikeInd_light(l)-padding:end, l, 1); NaN(abs(512-spikeInd_light(l)-padding),1)];
            lightSpikeWaveforms(:,l,2) = [spikeCSC_light(spikeInd_light(l)-padding:end, l, 2); NaN(abs(512-spikeInd_light(l)-padding),1)];
            lightSpikeWaveforms(:,l,3) = [spikeCSC_light(spikeInd_light(l)-padding:end, l, 3); NaN(abs(512-spikeInd_light(l)-padding),1)];
            lightSpikeWaveforms(:,l,4) = [spikeCSC_light(spikeInd_light(l)-padding:end, l, 4); NaN(abs(512-spikeInd_light(l)-padding),1)];
        end
    end
    for l = 1:length(spikeInd_spont)
        if spikeInd_spont(l) > padding && spikeInd_spont(l) < 512-padding
            spontSpikeWaveforms(:,l,1) = spikeCSC_spont(spikeInd_spont(l) - padding:spikeInd_spont(l) + padding, l, 1);
            spontSpikeWaveforms(:,l,2) = spikeCSC_spont(spikeInd_spont(l) - padding:spikeInd_spont(l) + padding, l, 2);
            spontSpikeWaveforms(:,l,3) = spikeCSC_spont(spikeInd_spont(l) - padding:spikeInd_spont(l) + padding, l, 3);
            spontSpikeWaveforms(:,l,4) = spikeCSC_spont(spikeInd_spont(l) - padding:spikeInd_spont(l) + padding, l, 4);
        elseif spikeInd_spont(l) <= padding
            spontSpikeWaveforms(:,l,1) = [NaN(padding-spikeInd_spont(l)+1,1); spikeCSC_spont(1:spikeInd_spont(l)+padding, l, 1)];
            spontSpikeWaveforms(:,l,2) = [NaN(padding-spikeInd_spont(l)+1,1); spikeCSC_spont(1:spikeInd_spont(l)+padding, l, 2)];
            spontSpikeWaveforms(:,l,3) = [NaN(padding-spikeInd_spont(l)+1,1); spikeCSC_spont(1:spikeInd_spont(l)+padding, l, 3)];
            spontSpikeWaveforms(:,l,4) = [NaN(padding-spikeInd_spont(l)+1,1); spikeCSC_spont(1:spikeInd_spont(l)+padding, l, 4)];
        elseif spikeInd_spont(l) >= 512-padding
            spontSpikeWaveforms(:,l,1) = [spikeCSC_spont(spikeInd_spont(l)-padding:end, l, 1); NaN(abs(512-spikeInd_spont(l)-padding),1)];
            spontSpikeWaveforms(:,l,2) = [spikeCSC_spont(spikeInd_spont(l)-padding:end, l, 2); NaN(abs(512-spikeInd_spont(l)-padding),1)];
            spontSpikeWaveforms(:,l,3) = [spikeCSC_spont(spikeInd_spont(l)-padding:end, l, 3); NaN(abs(512-spikeInd_spont(l)-padding),1)];
            spontSpikeWaveforms(:,l,4) = [spikeCSC_spont(spikeInd_spont(l)-padding:end, l, 4); NaN(abs(512-spikeInd_spont(l)-padding),1)];
        end
    end
    
    
    %% plot everything
    
    rasters = figure; subplot(4,3,[1:6]); hold on; title(strcat(spikeFields{i},' ;  ', filename),'Interpreter','none')
    xlabel('Time (ms)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','RelSpikeStartTime',-tB,'XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    for l = 1:2:length(laserRast{1})
        plotShaded([laserRast{1}(l)-tB laserRast{1}(l+1)-tB],[0 0; 1+numTrains 1+numTrains],'b');
    end
    
    subplot(4,3,[7]); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 25]); xlim([0 size(latencyToResponse,2)+1])
    plot(nanmean(latencyToResponse), 'b');
    plot(nanmean(latencyToResponseSham), 'k');
    legend('laser','control');
    
    subplot(4,3,10); hold on;
    xlabel('Pulse'); ylabel('Response Fraction'); ylim([-0.1 1.1]); xlim([0 size(numResponses,2)+1])
    plot(sum(numResponses)/size(numResponses,1), 'b');
    plot(sum(numResponsesSham)/size(numResponsesSham,1), 'k');

    
    xpad = 72;
    for currCSC = 1:4
        subplot(4,3,8); hold on
        ylabel('Amplitude (uV)')
        if ~isempty(lightSpikeWaveforms)
            plot([1:size(lightSpikeWaveforms,1)]/32-32*padding/1000 + (xpad*padding/1000)*(currCSC-1/2), nanmean(sCon*lightSpikeWaveforms(:,:,currCSC)'),'b');
        end
        xlim([(xpad*padding/1000)*(0) (xpad*padding/1000)*(4)])
        subplot(4,3,11); hold on
        xlabel('Time (ms)'); ylabel('Amplitude (uV)')
        plot([1:size(spontSpikeWaveforms,1)]/32-32*padding/1000 + (xpad*padding/1000)*(currCSC-1/2), nanmean(sCon*spontSpikeWaveforms(:,:,currCSC)'),'k');
        xlim([(xpad*padding/1000)*(0) (xpad*padding/1000)*(4)])

        subplot(4,3,[9 12]); hold on
        xlabel('Time (ms)'); ylabel('Amplitude (uV)')
        if ~isempty(lightSpikeWaveforms)
            errorfill([1:size(lightSpikeWaveforms,1)]/32-32*padding/1000 + (xpad*padding/1000)*(currCSC-1/2), nanmean(sCon*lightSpikeWaveforms(:,:,currCSC)'), nanstd(sCon*lightSpikeWaveforms(:,:,currCSC)'),'b');
        end
        errorfill([1:size(spontSpikeWaveforms,1)]/32-32*padding/1000 + (xpad*padding/1000)*(currCSC-1/2), nanmean(sCon*spontSpikeWaveforms(:,:,currCSC)'), nanstd(sCon*spontSpikeWaveforms(:,:,currCSC)'),'k');
        xlim([(xpad*padding/1000)*(0) (xpad*padding/1000)*(4)])
    end
        %% Save Figures
    set(rasters, 'Position', get(0,'Screensize'))
    saveFigurePDF(rasters,[saveDir '\optoID_' spikeFields{i}])
end