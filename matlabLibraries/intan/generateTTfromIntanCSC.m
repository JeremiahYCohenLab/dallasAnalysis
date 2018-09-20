function generateTTfromIntanCSC(sessionString, varargin)

[root, sep] = currComputer();
p = inputParser;
% default parameters if none given
p.addParameter('RecordingType', 'session');
p.addParameter('RemoveLick', false);
p.addParameter('MedianRemoval', false);
p.addParameter('HighPassCutoffInHz', 300);
p.addParameter('SamplingFreq', 20000);
p.addParameter('ThresholdFactor', 3);
p.addParameter('RefractorySamples', 20);
p.addParameter('AnalyzeSpecificChans', []);
p.addParameter('NumberOfChans', 32);


p.parse(varargin{:});

pd = parseSessionString_oM(sessionString, root, sep);
ephysPath = [pd.ephysPath p.Results.RecordingType sep];
dirToSave = [pd.ephysPath 'spikeSort_' p.Results.RecordingType sep];

ephysDir = dir(ephysPath);
if isempty(ephysDir)
    error('No ephys folder in %s', sessionString)
end
if ~isdir(dirToSave)
    mkdir(dirToSave);
end

RHDmask = contains({ephysDir.name},'rhd');
fprintf('Analyzing %s\n', pd.sessionFolder)
digIn = [];
sampToSave = [];
featToSave = [];
ts = 0;
loopInd = 1;

if isempty(p.Results.AnalyzeSpecificChans)
    chans = 1:p.Results.NumberOfChans;
else
    chans = p.Results.AnalyzeSpecificChans;
end

for currRHD_ind = find(RHDmask)
    %SAVE THIS STRUCT AND OPEN IF IT ALREADY EXISTS
    fprintf('File number: %d of %d \n', loopInd, sum(RHDmask));
    filename = ephysDir(currRHD_ind).name;
    [digIn, traces, ~] = readIntan(ephysPath, filename);
    
    if p.Results.MedianRemoval == true
        traces = traces - median(traces);
    end
        
    for chanInd = chans
        Wn = p.Results.HighPassCutoffInHz / (p.Results.SamplingFreq/2);
        [b, a] = butter(2, Wn, 'high');
        filtTrace = filtfilt(b, a, -traces(chanInd,:));

        tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds
        ts_interp = ts:tSamp:ts + tSamp*(length(filtTrace) - 1);

        % threshold and median method from Rey, Pedreira, Quiroga (2015)
        thresh = p.Results.ThresholdFactor*round(median(abs(filtTrace))/0.6745);
%        locs = sort([peakseek(filtTrace, p.Results.RefractorySamples, thresh) peakseek(-filtTrace, p.Results.RefractorySamples, thresh)]); % look for a peak, avoid 32 samples (1ms at 32kHz)
        locs = peakseek(filtTrace, p.Results.RefractorySamples, thresh);
        
        allLocs = unique(locs, 'stable');
        allLocs(allLocs > length(ts_interp) - p.Results.RefractorySamples) = []; % remove spikes within 1ms of the end of recording
        allLocs(allLocs < p.Results.RefractorySamples) = []; % remove spikes within 1ms of the beginning of recording
        
        % remove lick artifact from TT; remove the 1ms preceding every lick event
        if p.Results.RemoveLick == true
            lickInds = [];
            lickOn = find((digIn(5,1:end-1) == 0 &  digIn(5,2:end) > 0) | ...
                (digIn(6,1:end-1) == 0 &  digIn(6,2:end) > 0 ));
            for i = 1:length(lickOn)
                lickIndsTmp{i} = find(allLocs < lickOn(i) & allLocs > lickOn(i) - 50);
                if ~isempty(lickIndsTmp{i})
                    [~, lickIndsMax] = max(filtTrace(allLocs(lickIndsTmp{i})));
                    lickInds(i) = lickIndsTmp{i}(lickIndsMax);
                else
                    lickInds(i) = NaN;
                end
            end
            lickInds = lickInds(~isnan(lickInds));
            allLocs(lickInds) = [];
        end
        
        while any(diff(allLocs) < p.Results.RefractorySamples) % while there is an overlap in peaks within 1ms
            for i = 1:length(allLocs) - 1
                if allLocs(i + 1) < allLocs(i) + p.Results.RefractorySamples
                    % save each trace and find where the best peak is
                    tmp(1, :) = filtTrace(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                    [~, tmpi] = max(max(tmp));
                    allLocs(i) = allLocs(i) + tmpi - 1;
                    allLocs(i + 1) = allLocs(i);
                end
            end
            allLocs = unique(allLocs, 'stable');
        end

%        fprintf('\tTotal of %0.1d spikes thresholded.\n', length(allLocs));

        locMat = NaN(length(allLocs), 32);      %for taking 32 samples around peak, has to be 32 for spikesort3D
        sampBack = round(1/3*32) - 1;
        sampFor = round(2/3*32);
        while allLocs(end) > (length(filtTrace) - sampFor)  %gets rid of spikes to close to the end
            allLocs = allLocs(1:end-1);
            locMat = locMat(1:end-1, :);
        end
        while allLocs(1) < sampBack  %gets rid of spikes too close to the beginning
            allLocs = allLocs(2:end);
            locMat = locMat(2:end, :);
        end        
        for i = 1:length(allLocs)
            locMat(i, :) = allLocs(i) - sampBack:allLocs(i) + sampFor;
        end

        clear sampToSaveTemp featToSaveTemp tsToSaveTemp
        sampToSaveTemp(:,1,:) = filtTrace(locMat)'*100;          %multiply to scale for spikesort3D
        sampToSaveTemp(:,2,:) = zeros(size(sampToSaveTemp(:,1,:)));
        sampToSaveTemp(:,3,:) = zeros(size(sampToSaveTemp(:,1,:)));
        sampToSaveTemp(:,4,:) = zeros(size(sampToSaveTemp(:,1,:)));
        featToSaveTemp(1,:) = max(filtTrace(locMat)'*100);
        featToSaveTemp = [featToSaveTemp; zeros(3, size(featToSaveTemp,2))]; 
        featToSaveTemp(5,:) = min(filtTrace(locMat)'*100);
        featToSaveTemp = [featToSaveTemp; zeros(3, size(featToSaveTemp,2))]; 

        if loopInd == 1
            sampToSave{chanInd} = sampToSaveTemp;
            featToSave{chanInd} = featToSaveTemp;
            tsToSave{chanInd} = ts_interp(allLocs);        
        else
            sampToSave{chanInd} = cat(3, sampToSave{chanInd}, sampToSaveTemp);
            featToSave{chanInd} = [featToSave{chanInd}, featToSaveTemp];
            tsToSave{chanInd} = [tsToSave{chanInd} ts_interp(allLocs)];
        end
     
    end

    ts = ts_interp(end) + tSamp;
    loopInd = loopInd + 1;
    
end

for chanInd = chans
    fprintf('Saving file %d of %d \n', chanInd, chans(end))
    Mat2NlxSpike([dirToSave 'C_' num2str(chanInd) '.ntt'], 0, 1, [], [1 1 1 1 1], tsToSave{chanInd}, ...
        zeros(1, length(tsToSave{chanInd})), zeros(1, length(tsToSave{chanInd})), featToSave{chanInd}, sampToSave{chanInd});  
end
    
fprintf('Finished\n')