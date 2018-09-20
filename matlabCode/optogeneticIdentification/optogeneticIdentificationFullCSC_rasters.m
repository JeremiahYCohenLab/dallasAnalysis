function [latencyToResponse] = optogeneticIdentificationFullCSC_rasters(sessionsToImport, delayComp, TTdelay, CSCdelay)

animalName = strtok(sessionsToImport{:}, 'd'); animalName = animalName(2:end); 
sessionName = sessionsToImport{:};
if exist(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Sorted\laserStruct.mat'),'file')
    load(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Sorted\laserStruct.mat'))
else
    [~, ~, laserStruct] = generateTrialStruct_uU(sessionsToImport);
    save(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Sorted\laserStruct.mat'), 'laserStruct')
end

if delayComp == 0 % this compensation is for DSP delay compensation disabled sessions
    TTdelay = 0;
    CSCdelay = 0;
    spikeDelay = 0;
elseif delayComp == 1
    spikeDelay = TTdelay - CSCdelay;
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

for i = clust'
    
    TTchann = str2double(spikeFields{i}(3));
    CSCchann{1} = ['CSC' num2str(4*TTchann-3)];
    CSCchann{2} = ['CSC' num2str(4*TTchann-2)]; 
    CSCchann{3} = ['CSC' num2str(4*TTchann-1)];
    CSCchann{4} = ['CSC' num2str(4*TTchann-0)];
    trainInd = 1:numPulses:length(laserStruct.laserOn);
    
    latencyToResponse = []; % latency of light-evoked responses
    numResponses = []; % success of light-evoked response
    
    spikeTimes_light_TT = []; % light-evoked spike in TT time
    spikeTimes_spont_TT = [];
    
    allCSCtimes = [];
    allCSCdata = [];

    
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
                latencyToResponse(ceil(j/numPulses),k) = tempLat - TTdelay/1000;
                numResponses(ceil(j/numPulses),k) = 1;
                spikeTimes_light_TT = [spikeTimes_light_TT tempSpike];
            else
                latencyToResponse(ceil(j/numPulses),k) = NaN;
                numResponses(ceil(j/numPulses),k) = 0;
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

    [allCSCtimes, allCSCdata(:,:,1)] = Nlx2MatCSC(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Unsorted\', CSCchann{1},'.ncs'),[1 0 0 0 1], 0, 1);
    [allCSCdata(:,:,2)] = Nlx2MatCSC(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Unsorted\', CSCchann{2},'.ncs'),[0 0 0 0 1], 0, 1);
    [allCSCdata(:,:,3)] = Nlx2MatCSC(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Unsorted\', CSCchann{3},'.ncs'),[0 0 0 0 1], 0, 1);
    [allCSCdata(:,:,4)] = Nlx2MatCSC(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, '\Unsorted\', CSCchann{4},'.ncs'),[0 0 0 0 1], 0, 1);
    allCSCdata = reshape(allCSCdata, [], 1, 4); allCSCdata = squeeze(allCSCdata);
    
    spikeInd_light = round(((1000*spikeTimes_light_TT - spikeDelay) - allCSCtimes(1))/31.25);
    spikeInd_spont = round(((1000*spikeTimes_spont_TT - spikeDelay) - allCSCtimes(1))/31.25);
    allCSCtimes = (allCSCtimes(1):31.25:allCSCtimes(end)) - allCSCtimes(1);

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
    
    figure; subplot(4,3,[1:6]); hold on; title(strcat(spikeFields{i},' ;  ', sessionName),'Interpreter','none')
    xlabel('Time (ms)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','RelSpikeStartTime',-tB,'XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    for l = 1:2:length(laserRast{1})
        plotShaded([laserRast{1}(l)-tB laserRast{1}(l+1)-tB],[0 0; 1+numTrains 1+numTrains],'b');
    end
    
    subplot(4,3,[7]); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 25]); xlim([0 size(latencyToResponse,2)+1])
    plot(nanmean(latencyToResponse));
    
    subplot(4,3,10); hold on;
    xlabel('Pulse'); ylabel('Response Fraction'); ylim([-0.1 1.1]); xlim([0 size(numResponses,2)+1])
    plot(sum(numResponses)/size(numResponses,1));
    
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
end