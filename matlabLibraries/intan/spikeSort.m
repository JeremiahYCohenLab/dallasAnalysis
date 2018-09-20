% Spike sort data acquired from Intan system
% 'ampFiles' is cell array of names of tetrode channels
% Written by Jeremiah Y. Cohen
%
% TODO: configure for stereotrode or single-electrode configurations

function spikeTimes = spikeSort(ampFiles, path, samplingRate,...
    triggeredWindowSamples, ISIviol)
    
    filtData = cell(0);
    thresh = [];
    [b, a] = butter(4, [.01 .9], 'bandpass');
    
    % Raw signals
    for i = 1 : 4
        fileinfo = dir(strcat(path, ampFiles{i}));
        numSamples = fileinfo.bytes/2; % int16 = 2 bytes
        fid = fopen(strcat(path, ampFiles{i}), 'r');
        rawSignal = fread(fid, numSamples, 'int16');
        fclose(fid);
        % convert to microvolts
        rawSignal = rawSignal * 0.195;
        
        % High-pass (zero-phase) filter data
%        dataIn = rawSignals{i};
        filtData = [filtData; filtfilt(b, a, rawSignal)];
        
        % Plot a subset of each channel for choosing thresholds
        figure; hold on;
        plot(0:(1/samplingRate):10, rawSignal(1:(1+samplingRate*10)));
        plot(0:(1/samplingRate):10, ...
            filtData{i}(1:(1+samplingRate*10)), 'r');
        [~, y] = getpts;
        line([0 10], [y y]);
        hold off;
        thresh = [thresh; y];
        %replotTF = uicontrol('Style', 'checkbox', 'String', 'Done');
        
        % if threshold is negative, reflect signal across horizontal axis
        if (thresh(i) < 0)
            filtData{i} = abs(filtData{i} - max(...
                filtData{i})) - max(filtData{i});
        end
    end
    
    spikeTimes = cell(0);
    whichSpikeTimes = cell(0);
    peakTimes = [];
    triggeredWaveforms = cell(0);
    triggeredPeaks = cell(0);
    triggeredValleys = cell(0);
    %peakVals = cell(0);
    %tempTimes = [];
    for i = 1 : 4
        [peakVals, temp] = findpeaks(filtData{i}, 'MinPeakHeight',...
            abs(thresh(i)), 'MinPeakDistance', ISIviol, 'Threshold', 1e-4);
        
        % remove peaks earlier than half triggered window size from start
        peakVals = peakVals(temp > triggeredWindowSamples/3);
        temp = temp(temp > triggeredWindowSamples/3);
        
        % remove peaks later than half triggered window size from end
        peakVals = peakVals(temp < length(...
            filtData{i}) - 2*triggeredWindowSamples/3);
        peakTimes = [peakTimes; temp(temp < length(...
            filtData{i}) - 2*triggeredWindowSamples/3)];

        %tempTimes = [tempTimes; temp];
    end
    
    peakTimes = unique(sort(peakTimes));
    
    for i = 1 : 4
        % Plot all triggered waveforms
        triggeredWaveforms = [triggeredWaveforms; ...
            zeros(length(peakTimes), triggeredWindowSamples)];
        for j = 1 : length(peakTimes)
            triggeredWaveforms{i}(j, :) = filtData{i}(...
                (1 + peakTimes(j) - triggeredWindowSamples/3) : (...
                peakTimes(j) + 2*triggeredWindowSamples/3));
        end
        
        figure; hold on;
        plot(triggeredWaveforms{i}(1, :));
        for j = 2 : length(peakTimes)
            plot(triggeredWaveforms{i}(j, :));
        end
        hold off;
        % calculate peak and valley voltages
        triggeredPeaks = [triggeredPeaks; ...
            max(triggeredWaveforms{i}, [], 2)];
        triggeredValleys = [triggeredValleys; ...
            min(triggeredWaveforms{i}, [], 2)];
    end
    
    for i = 1 : 4
        temp = triggeredPeaks{i};
        triggeredPeaks{i}(temp > 2000) = 0;
    end
    
    figure; hold on; iter = 0;
    for i = 1 : 4
        for j = 1 : 4
            iter = iter + 1;
            subplot(4, 4, iter);
            hexscatter(triggeredPeaks{i}, triggeredPeaks{j}, 'res', 200)
        end
    end
    hold off;
    figure; hold on; iter = 0;
    for i = 1 : 4
        for j = 1 : 4
            iter = iter + 1;
            subplot(4, 4, iter);
            hexscatter(triggeredValleys{i}, triggeredValleys{j}, 'res',200)
        end
    end
    hold off;
    figure; hold on; iter = 0;
    for i = 1 : 4
        for j = 1 : 4
            iter = iter + 1;
            subplot(4, 4, iter);
            hexscatter(triggeredPeaks{i}, triggeredValleys{j}, 'res', 200)
        end
    end
    hold off;
    
    sortTF = true;
    while sortTF
        % done sorting?
        sortTF = str2num(cell2mat(inputdlg(...
            {'Still sorting this tetrode?'}, 'Input', 1, {'false'})));
        
        % choose features on which to cluster
        prompt = {'Enter first feature (1 = peak1, 8 = valley4):'; ...
            'Enter second feature:'};
        features = str2num(cell2mat(inputdlg(prompt,'Input',1,{'1','2'})));

        if features(1) < 5
            myX = triggeredPeaks{features(1)};
        else
            myX = triggeredValleys{features(1)-4};
        end
        if features(2) < 5
            myY = triggeredPeaks{features(2)};
        else
            myY = triggeredValleys{features(2)-4};
        end
    
        % plot features and select cluster
        figure; hold on;
        hexscatter(myX, myY, 'res', 200)
        [x, y] = getpts;
        insideRegion = inpolygon(myX, myY, x, y);
        plot(myX(insideRegion), myY(insideRegion), 'r.')
        hold off
        
        spikeTimes = [spikeTimes; {peakTimes(insideRegion)}];
        whichSpikeTimes = [whichSpikeTimes; {find(insideRegion)}];
        
        % plot cluster on all feature pairs
        figure; hold on; iter = 0;
        for i = 1 : 4
            for j = 1 : 4
                iter = iter + 1;
                subplot(4, 4, iter); hold on;
                hexscatter(triggeredPeaks{i},triggeredPeaks{j},'res',200)
                plot(triggeredPeaks{i}(insideRegion), ...
                    triggeredPeaks{j}(insideRegion), 'r.')
                hold off;
            end
        end
        hold off;
        figure; hold on; iter = 0;
        for i = 1 : 4
            for j = 1 : 4
                iter = iter + 1;
                subplot(4, 4, iter); hold on;
                hexscatter(triggeredValleys{i},...
                    triggeredValleys{j},'res',200)
                plot(triggeredValleys{i}(insideRegion), ...
                    triggeredValleys{j}(insideRegion), 'r.')
                hold off;
            end
        end
        hold off;
        figure; hold on; iter = 0;
        for i = 1 : 4
            for j = 1 : 4
                iter = iter + 1;
                subplot(4, 4, iter); hold on;
                hexscatter(triggeredPeaks{i},triggeredValleys{j},'res',200)
                plot(triggeredPeaks{i}(insideRegion), ...
                    triggeredValleys{j}(insideRegion), 'r.')
                hold off;
            end
        end
        hold off;

    end
    
    % plot all waveforms
    for i = 1 : length(spikeTimes)
        figure;
        for j = 1 : 4
            subplot(2, 2, j); hold on;
            plot(triggeredWaveforms{j}(whichSpikeTimes{i}(1), :));
            for k = 2 : length(whichSpikeTimes{i})
                plot(triggeredWaveforms{j}(whichSpikeTimes{i}(k), :));
            end
            hold off;
        end
    end