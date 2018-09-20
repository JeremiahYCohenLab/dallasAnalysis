function Nlx_recenterTTfromCSCtemp(sessionString, varargin)

[root, sep] = currComputer_operantMatching();
p = inputParser;
% default parameters if none given
p.addParameter('Root', root);
p.addParameter('Separator', sep)
p.addParameter('SaveFolder', 'C:\Users\Cooper\Desktop\');
p.addParameter('RecordingType', 'session')
p.addParameter('RemoveLick_Flag', false)
p.addParameter('HighPassCutoffInHz', 300)
p.addParameter('SamplingFreq', 32000)
p.addParameter('ThresholdFactor', 4);
p.addParameter('RefractorySamples', 32);
p.addParameter('AnalyzeSpecificTTs', []);
% p.addParameter('CSCscaleFactor', []);
% p.addParameter('CSCstoScale', []);

p.parse(varargin{:});

pd = parseSessionString_oM(sessionString, p.Results.Root, p.Results.Separator);
if regexp(p.Results.RecordingType, 'session') & isstrprop(pd.suptitleName(end), 'alpha')
    ephysPath = [pd.ephysPath 'session ' pd.suptitleName(end) sep];
elseif regexp(p.Results.RecordingType, 'opto')
    ephysPath = [pd.ehysPath 'opto' sep];
else
    ephysPath = [pd.ephysPath 'session' sep];
end
ephysDir = dir(ephysPath);
if isempty(ephysDir)
    error('No ephys folder in %s', sessionString)
end
TTmask = contains({ephysDir.name},'TT');

dirToSave = [pd.ephysPath p.Results.Separator 'spikeSort' p.Results.Separator];
if ~isdir(dirToSave)
    mkdir(dirToSave);
end

% % test for pauses in session
% ts = Nlx2MatCSC([ephysPath 'CSC1.ncs'], [1 0 0 0 0], 0, 1, []);
% tsDifferences = unique(diff(ts), 'stable');
% tsDifferences = tsDifferences(2:end);
% if ~isempty(tsDifferences)
%     fprintf('Session was paused %i times.\n', length(tsDifferences))
%     CSCstoScale = input('Which CSCs should be scaled? Enter nothing to apply scaling uniformly. ');
%     CSCscaleFactor = input('What should the scale factors be? Enter nothing to scale nothing. ');
%     if ~isempty(CSCscaleFactor)
%         if length(tsDifferences) ~= length(CSCscaleFactor)
%             error('length(tsDifferences) (%i) and length(CSCscaleFactor) (%i) do not match.\n\n', length(tsDifferences), length(CSCscaleFactor))
%         else % lengths match
%             if isempty(CSCstoScale)
%                 fprintf('Session paused but no specific CSCs given; assuming all will be scaled by CSCscaleFactor %s.\n\n', num2str(CSCscaleFactor))
%             else
%                 fprintf('Only scaling the following CSCs: %s\n\n', num2str(CSCstoScale))
%             end
%         end
%     else
%         fprintf('Session was paused but no CSCscaleFactor inputted; assuming no changes in scaling.\n\n')
%     end
% else
%     fprintf('No pausing in session; assuming no rescaling necessary.\n\n')
% end

fprintf('Analyzing %s\n', pd.sessionFolder)
for currTT_ind = find(TTmask)
    TTname = ephysDir(currTT_ind).name;
    TTnum = TTname(3:regexp(TTname,'.ntt') - 1);
    TTnum = str2double(TTnum);
    
    if isempty(p.Results.AnalyzeSpecificTTs) || (~isempty(p.Results.AnalyzeSpecificTTs) && any(TTnum == p.Results.AnalyzeSpecificTTs))
        chan0 = TTnum + (3*(TTnum - 1));
        chan1 = chan0 + 1;
        chan2 = chan0 + 2;
        chan3 = chan0 + 3;

        fprintf('Currently on TT%0.1d: CSC %0.1d,%0.1d,%0.1d,%0.1d. ', TTnum, chan0, chan1, chan2, chan3)
        [ts, samp0, head0] = Nlx2MatCSC([ephysPath 'CSC' num2str(chan0) '.ncs'], [1 0 0 0 1], 1, 1, []);
        samp1 = Nlx2MatCSC([ephysPath 'CSC' num2str(chan1) '.ncs'], [0 0 0 0 1], 0, 1, []);
        samp2 = Nlx2MatCSC([ephysPath 'CSC' num2str(chan2) '.ncs'], [0 0 0 0 1], 0, 1, []);
        samp3 = Nlx2MatCSC([ephysPath 'CSC' num2str(chan3) '.ncs'], [0 0 0 0 1], 0, 1, []);
        
%         if ~isempty(tsDifferences) % any pauses in this session; indicative of change in CSC
%             indMin = 1;
%             for currPause = 1:length(tsDifferences)
%                 indMax = find(diff(ts) == tsDifferences(currPause));
%                 scaleFactor = CSCscaleFactor(currPause);
%                 if isempty(CSCstoScale) % rescale all CSCs
%                     samp0(:, indMin:indMax) = scaleFactor*samp0(:, indMin:indMax);
%                     samp1(:, indMin:indMax) = scaleFactor*samp1(:, indMin:indMax);
%                     samp2(:, indMin:indMax) = scaleFactor*samp2(:, indMin:indMax);
%                     samp3(:, indMin:indMax) = scaleFactor*samp3(:, indMin:indMax);
%                     fprintf('\t Rescaling CSCs%i,%i,%i,%i from indices %i to %i\n', chan0, chan1, chan2, chan3, indMin, indMax)
%                 else
%                     if any(CSCstoScale == chan0)
%                         samp0(:, indMin:indMax) = scaleFactor*samp0(:, indMin:indMax);
%                         fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan0, indMin, indMax)
%                     end
%                     if any(CSCstoScale == chan1)
%                         samp1(:, indMin:indMax) = scaleFactor*samp1(:, indMin:indMax);
%                         fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan1, indMin, indMax)
%                     end
%                     if any(CSCstoScale == chan2)
%                         samp2(:, indMin:indMax) = scaleFactor*samp2(:, indMin:indMax);
%                         fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan2, indMin, indMax)
%                     end
%                     if any(CSCstoScale == chan3)
%                         samp3(:, indMin:indMax) = scaleFactor*samp3(:, indMin:indMax);
%                         fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan3, indMin, indMax)
%                     end
%                 end
%                 indMin = indMax + 1;
%             end
%         end



        Wn = p.Results.HighPassCutoffInHz / (p.Results.SamplingFreq/2);
        [b, a] = butter(2, Wn, 'high');
        samp0 = filtfilt(b, a, samp0(:));
        samp1 = filtfilt(b, a, samp1(:));
        samp2 = filtfilt(b, a, samp2(:));
        samp3 = filtfilt(b, a, samp3(:));

        tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(length(samp0) - 1);

        % threshold and median method from Rey, Pedreira, Quiroga (2015)
        thresh0 = p.Results.ThresholdFactor*round(median(abs(samp0))/0.6745);
        thresh1 = p.Results.ThresholdFactor*round(median(abs(samp1))/0.6745);
        thresh2 = p.Results.ThresholdFactor*round(median(abs(samp2))/0.6745);
        thresh3 = p.Results.ThresholdFactor*round(median(abs(samp3))/0.6745);
        locs0 = peakseek(samp0, p.Results.RefractorySamples, thresh0); % look for a peak, avoid 32 samples (1ms at 32kHz)
        locs1 = peakseek(samp1, p.Results.RefractorySamples, thresh1);
        locs2 = peakseek(samp2, p.Results.RefractorySamples, thresh2);
        locs3 = peakseek(samp3, p.Results.RefractorySamples, thresh3);

        allLocs = unique(sort([locs0 locs1 locs2 locs3]), 'stable');
        allLocs(allLocs > length(ts_interp) - p.Results.RefractorySamples) = []; % remove spikes within 1ms of the end of recording
        allLocs(allLocs < p.Results.RefractorySamples) = []; % remove spikes within 1ms of the beginning of recording
        while any(diff(allLocs) < p.Results.RefractorySamples) % while there is an overlap in peaks within 1ms
            for i = 1:length(allLocs) - 1
                if allLocs(i + 1) < allLocs(i) + p.Results.RefractorySamples
                    % save each trace and find where the best peak is
                    tmp(1, :) = samp0(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                    tmp(2, :) = samp1(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                    tmp(3, :) = samp2(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                    tmp(4, :) = samp3(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                    [~, tmpi] = max(max(tmp));
                    allLocs(i) = allLocs(i) + tmpi - 1;
                    allLocs(i + 1) = allLocs(i);
                end
            end
            allLocs = unique(allLocs, 'stable');
        end

        % remove lick artifact from TT; remove the 1ms preceding every lick event
        if p.Results.RemoveLick_Flag == true
            [tsTTLs, TTLs] = Nlx2MatEV([ephysPath 'Events.nev'], [1 0 1 0 0], 0, 1);
            lickTTLs = [16 32 18 33 82 97]; % possible lick events on 295H Neuralynx rig
            tsLicks = tsTTLs(ismember(TTLs, lickTTLs));
            lickInds = find(ismember(round(tsLicks/10)*10, round(ts_interp/10)*10)); % round due to slide differences in times from .nev file and .ncs file

            for i = 1:length(lickInds)
                tmpLickArtifact = allLocs - lickInds(i);
                allLocs(tmpLickArtifact >= -p.Results.SamplingFreq*1e-3 & tmpLickArtifact < 0) = [];
            end
        end

        fprintf('\tTotal of %0.1d spikes thresholded.\n', length(allLocs));

        locMat = NaN(length(allLocs), p.Results.RefractorySamples);
        sampBack = round(1/3*32) - 1;
        sampFor = round(2/3*32);
        for i = 1:length(allLocs)
            locMat(i, :) = allLocs(i) - sampBack:allLocs(i) + sampFor;
        end

        clear sampToSave featToSave tsToSave
        sampToSave(:,1,:) = samp0(locMat)';
        sampToSave(:,2,:) = samp1(locMat)';
        sampToSave(:,3,:) = samp2(locMat)';
        sampToSave(:,4,:) = samp3(locMat)';

        featToSave(1,:) = max(samp0(locMat)');
        featToSave(2,:) = max(samp1(locMat)');
        featToSave(3,:) = max(samp2(locMat)');
        featToSave(4,:) = max(samp3(locMat)');
        featToSave(5,:) = min(samp0(locMat)');
        featToSave(6,:) = min(samp1(locMat)');
        featToSave(7,:) = min(samp2(locMat)');
        featToSave(8,:) = min(samp3(locMat)');

        tsToSave = ts_interp(allLocs);

        Mat2NlxSpike([dirToSave TTname], 0, 1, [], [1 1 1 1 1], tsToSave, ...
            zeros(1, length(tsToSave)), zeros(1, length(tsToSave)), featToSave, sampToSave);
    end
end

fprintf('Finished\n')