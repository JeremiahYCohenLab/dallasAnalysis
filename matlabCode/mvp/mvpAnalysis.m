function [traceTbl, spikeStruct, rawSpikeStruct, maxCrossSpikes] = mvpAnalysis(xlFile, sheet)

filesLoc = 'C:\Users\cooper_PC\Desktop\mvp\Eastman Matlab for Cooper\Cells\';

[numbers, dates, ~] = xlsread(xlFile, sheet);
cellNumbers = numbers(:,1);
fourAP = numbers(:,2);
hyperpolCurrent = -numbers(:,3);
dates = dates(2:end,1);
numCells = length(dates);
cellNames = cell(1, numCells);

%set variables
plotFile = 1;
latStartThresh = -70;
sampRate = 20000;
hyperpolWindow = [20001:25000];
hyperpolTauWindow = [(hyperpolWindow(1) + 0.1*length(hyperpolWindow)) : (hyperpolWindow(end) - 0.1*length(hyperpolWindow))];
depolWindow = [25101:30000];
maxSpikes = 10;

spikeStruct = struct;                   %make structures for trace features across cells and sweeps
rawSpikeStruct = struct;

restVm = cell(numCells,1);          %make matrices for trace features
hyperpolTau = cell(numCells,1);         
hyperpolMax = cell(numCells,1);
ir = cell(numCells,1);
firstSpikeLat = cell(numCells,1);
spikeFreq = cell(numCells,1);
spikeNum = cell(numCells,1);
firstSpikeHW = cell(numCells,1);

maxCrossSpikes = [];

for cellInd = 1:numCells
    cellName = sprintf([dates{cellInd} ' Cell %i'], cellNumbers(cellInd));
    cellPath = [filesLoc cellName '\'];
    files = dir(cellPath);
    if fourAP(cellInd) == 1
        sweepFileInds = find(contains({files.name},'.txt') & contains({files.name},'4ap'));
        cellName = [cellName '_4ap'];
    else
        sweepFileInds = find(contains({files.name},'.txt') & ~contains({files.name},'4ap'));
    end
    sweepFileNames = {files(sweepFileInds).name};
    numSweeps = length(sweepFileNames);
    
    spikeThresh = nan(maxSpikes,numSweeps);        %make matrices for spike features of variable number of spikes
    spikeAmp = nan(maxSpikes,numSweeps);
    spikePeak = nan(maxSpikes,numSweeps);
    spikeHW = nan(maxSpikes,numSweeps);            %half width
    for i = 1:maxSpikes
        rowNames{i} = sprintf(['spike_%i'], i);
    end
        
 
    hyperpolTau{cellInd} = nan(1,numSweeps);         %make matrices for trace features
    hyperpolMax{cellInd} = nan(1,numSweeps);
    ir{cellInd} = nan(1,numSweeps);
    firstSpikeLat{cellInd} = nan(1,numSweeps);
    spikeFreq{cellInd} = nan(1,numSweeps);
    spikeNum{cellInd} = nan(1,numSweeps);
    firstSpikeHW{cellInd} = nan(1,numSweeps);
    
    %run analysis on each sweep file
    for i = 1:numSweeps
        fid = fopen([cellPath sweepFileNames{i}]);
        C(1,i) = textscan(fid, '%f');
        trace = C{1,i};

        allSpikeTimes = peakseek(trace, 100, 0);         %find peaks with min distance 100 samples (5ms) and min peak 0 mV
        spikeTimes = allSpikeTimes(find(allSpikeTimes > hyperpolWindow(1)+500 & allSpikeTimes < depolWindow(end)));
        spikeNumInd = length(spikeTimes);
        spikeNum{cellInd}(i) =  spikeNumInd;
        
        %find second derivative to find spike beginning
        potentialChangeRate = gradient(gradient(trace));                                % find rate of change of potential
        potentialChangeRate_threshCross = find(potentialChangeRate(depolWindow) > 0.1);     % find indices that exceed slope threshold of 40mV/s^2
        spikeThreshTimes = potentialChangeRate_threshCross + depolWindow(1) - 1;                   % find spike starts in overall time
        
        %delete multiple instances of threshold crossing on the same spike
        spikeThreshInds = [1 NaN(1,length(spikeTimes)-1)];
        for j = 1:length(spikeTimes)-1
           indTmp = find(spikeThreshTimes > spikeTimes(j+1) - 100 & spikeThreshTimes < spikeTimes(j+1), 1);     %find thresh crossings within 5ms before spike peak
           if isempty(indTmp) %if there aren't any threshold crossings, find max 2nd deriv crossing
               [potentialChangeRateMax,spikeThreshTimeTmp] = max(potentialChangeRate(spikeTimes(j+1)-100:spikeTimes(j+1)));
               spikeThreshTimeTmp = spikeThreshTimeTmp + spikeTimes(j) - 1;
               spikeThreshTimes = sort([spikeThreshTimes; spikeThreshTimeTmp]);
               indTmp = find(spikeThreshTimes == spikeThreshTimeTmp);
               maxCrossSpikes = [maxCrossSpikes; {sprintf([cellName ' sweep %i spike %i'], [i j])}...
                   (potentialChangeRateMax / (1/sampRate)^2) / (1000^2)]; %make table of spikes where this occurs, and max rate of change
           end
           spikeThreshInds(j+1) = indTmp;
        end
            
        spikeThreshTimes = spikeThreshTimes(spikeThreshInds);
        spikeThresh(1:spikeNumInd,i) = trace(spikeThreshTimes);
        
        spikePeak(1:spikeNumInd,i) = trace(spikeTimes);                                    %find peak potential of each spike
        spikeAmp(1:spikeNumInd,i) = spikePeak(1:spikeNumInd,i) - spikeThresh(1:spikeNumInd,i);   %find peak amplitude
        
        %find spike frequency between first two spikes
        if spikeNumInd > 1
            isiTmp = diff(spikeTimes) / sampRate;
            spikeFreq{cellInd}(i) = 1 / isiTmp(1);
        end
        
        %find half-width
        for j = 1:spikeNumInd
            halfAmp = (spikeAmp(j,i) / 2) + spikeThresh(j,i);                            %find half amplitude
            halfRiseInd = find(trace(spikeThreshTimes(j):spikeTimes(j)) > halfAmp, 1);   %find point where trace crosses it
            x = [halfRiseInd-1 halfRiseInd];                                             %find line between points on sides of threshold
            x = x + spikeThreshTimes(j) - 1;
            y = trace(x);
            c = [[1; 1]  x(:)]\y(:);                                                        %calculate parameter vector
            halfRise = (halfAmp - c(1)) / c(2);                                             %x = (y-b)/ m
            halfFallInd = find(trace(spikeTimes(j):end) < halfAmp, 1);
            x = [halfFallInd-1 halfFallInd];
            x = x + spikeTimes(j) - 1;
            y = trace(x);
            c = [[1; 1]  x(:)]\y(:);            
            halfFall = (halfAmp - c(1)) / c(2);
            spikeHW(j,i) = (halfFall - halfRise) / sampRate * 1000;                      %find difference and convert to time in ms
        end
        
        traceInd = 101;
        while any(abs(traceInd - allSpikeTimes) < 100)
            traceInd = traceInd + 100;
        end
        if traceInd > hyperpolWindow(1)
            restVm{cellInd}(i) = NaN;
        else
            restVm{cellInd}(i) = mean(trace(traceInd-100:traceInd+100));
        end
            
        
        hyperpolMax{cellInd}(i) = min(trace(hyperpolWindow));            %find max hyperpolarization potential
        ir{cellInd}(i) = (hyperpolMax{cellInd}(i) / hyperpolCurrent(cellInd)) * 1000;      %find input resistance in MOhms
        
        %find tau of hyperpolarization
        hyperpolAmp = hyperpolMax{cellInd}(i) - trace(hyperpolWindow(1));
        hyperpolInds(1) = find(trace(hyperpolWindow) < (trace(hyperpolWindow(1)) + 0.1 * hyperpolAmp), 1);
        hyperpolInds(2) = find(trace(hyperpolWindow) < (trace(hyperpolWindow(1)) + 0.9 * hyperpolAmp), 1);
        x = [1:diff(hyperpolInds)+1] / sampRate' * 1000;
        ft = fittype( 'a*exp((-1/b)*x) - c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares'); opts.Display = 'Off'; opts.StartPoint = [50 20 75];
        hyperpolFit = fit(x', trace(hyperpolWindow(hyperpolInds(1):hyperpolInds(2))), ft, opts );
        fitCoeffs = coeffvalues(hyperpolFit); 
        hyperpolTau{cellInd}(i) = fitCoeffs(2);

        
        %find starting point for latency to first spike calculation
        latStartInd = find(trace(hyperpolWindow(end):spikeThreshTimes(1)) > latStartThresh, 1);
        x = [latStartInd-1 latStartInd];                                             %find line between points on sides of threshold
        x = x + hyperpolWindow(end) - 1;
        y = trace(x);
        c = [[1; 1]  x(:)]\y(:);                                                        %calculate parameter vector
        latStart = (latStartThresh - c(1)) / c(2);                                             %x = (y-b)/ m
        firstSpikeLat{cellInd}(i) = (spikeThreshTimes(1) - latStart) / sampRate * 1000;             %find difference and convert to time in ms
        
        fclose(fid);    %close open text file
    end
    
    %put data for every spike and sweep in a structure
    cellName = erase(cellName, '-'); cellName = erase(cellName, ' '); cellName = ['d' cellName];
    cellNames{cellInd} = cellName;
    
    rawSpikeStruct.(cellName).spikeThresh = table(spikeThresh, 'VariableNames', {'sweeps'}, 'RowNames', rowNames);
    rawSpikeStruct.(cellName).spikeAmp = table(spikeAmp, 'VariableNames', {'sweeps'}, 'RowNames', rowNames);
    rawSpikeStruct.(cellName).spikePeak = table(spikePeak, 'VariableNames', {'sweeps'}, 'RowNames', rowNames);
    rawSpikeStruct.(cellName).spikeHW = table(spikeHW, 'VariableNames', {'sweeps'}, 'RowNames', rowNames);
    
    %put average and sd of data for each spike in a structure
    rowNames = {};
    for i = 1:length(spikeThresh)
        rowNames{i} = sprintf('spike_%i',i);
    end
    tblTmp = table([nanmean(spikeThresh,2) nanstd(spikeThresh,0,2)], [nanmean(spikeAmp,2) nanstd(spikeAmp,0,2)],...
        [nanmean(spikePeak,2) nanstd(spikePeak,0,2)], [nanmean(spikeHW,2) nanstd(spikeHW,0,2)], 'VariableNames',...
        {'spikeThresh' 'spikeAmp' 'spikePeak' 'spikeHW'}, 'RowNames', rowNames);
    spikeStruct.(cellName) = tblTmp;
    
    
    %find average and sd of trace features
    restVm{cellInd} = [mean(restVm{cellInd}) nanstd(restVm{cellInd},0)];
    hyperpolTau{cellInd} = [mean(hyperpolTau{cellInd}) nanstd(hyperpolTau{cellInd},0)];
    hyperpolMax{cellInd} = [mean(hyperpolMax{cellInd}) nanstd(hyperpolMax{cellInd},0)];
    ir{cellInd} = [mean(ir{cellInd}) nanstd(ir{cellInd},0)];
    firstSpikeLat{cellInd} = [mean(firstSpikeLat{cellInd}) nanstd(firstSpikeLat{cellInd},0)];
    spikeFreq{cellInd} = [mean(spikeFreq{cellInd}) nanstd(spikeFreq{cellInd},0)];
    spikeNum{cellInd} = [mean(spikeNum{cellInd}) nanstd(spikeNum{cellInd},0)];
    firstSpikeHW{cellInd} = [nanmean(spikeHW(1,:),2) nanstd(spikeHW(1,:),0,2)];
    
end

traceTbl = table(cell2mat(restVm), cell2mat(hyperpolTau), cell2mat(hyperpolMax), cell2mat(ir), cell2mat(firstSpikeLat),...
    cell2mat(spikeFreq), cell2mat(spikeNum), cell2mat(firstSpikeHW), 'VariableNames',...
    {'restVm' 'hyperpolTau' 'hyperpolMax' 'ir' 'firstSpikeLat' 'spikeFreq' 'spikeNum' 'firstSpikeHW'}, 'RowNames', cellNames);


end