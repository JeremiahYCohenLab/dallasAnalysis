function filterCSCtraces(tracePath, CSCchann, spikeTimes_light_TT)

if ~exist(strcat(tracePath, 'filtered\'), 'dir')
    mkdir(strcat(tracePath, 'filtered\'));
end

[header] = Nlx2MatCSC(strcat(tracePath, CSCchann{1}, '.ncs'), [0 0 0 0 0], 1, 1);       %get header info from trace
[cell,~] = find(~cellfun(@isempty,strfind(header, '-SamplingFrequency')) == 1);
[~, fs] = strtok(header(cell), ' ');                                                    %retrieve sampling frequency

fs = str2double(fs{1}(2:end));  %sample rate in Hz
fcutlow = 15999;   %low cut frequency in Hz
order = 16;   %order of filter
[b,a] = butter(order,fcutlow/(fs/2),'high');   %butterworth bandpass filter design

recordSize = 512;

for i = 1:length(CSCchann)

    [timeStamps,trace, header] = Nlx2MatCSC(strcat(tracePath, CSCchann{i}, '.ncs'), [1 0 0 0 1], 1, 1);
    filtTrace = reshape(trace', [1,(size(trace, 2)*recordSize)]);
    filtTrace = filtfilt(b,a,filtTrace);
    outTrace = reshape(filtTrace, [length(filtTrace)/recordSize recordSize])';
    
    Mat2NlxCSC(strcat(tracePath, 'filtered\', CSCchann{i}, '.ncs'), 0, 1, [], [1 0 0 0 1 1], timeStamps, outTrace, header);
    
end
