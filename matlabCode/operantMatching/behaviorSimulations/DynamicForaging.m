classdef DynamicForaging
    %
    %   DynamicForaging class
    %
    % A class that simulates the VI-VI dynamic foraging task
    % Initialization: Takes as inputs reward probabilities, block lengths, and number of trials
    % Running: Takes as inputs the choice
    %
    % BAB 6/21/17
    %
    properties
        BaitLeft_Flag = false;
        BaitRight_Flag = false;
        RewardProbabilities % left, right
        BlockLength
        BlockEnd
        AllRewards = [];
        AllChoices = [];
        NewBlock_Flag = false;
        BlockSwitch = [];
        Trial
        MaxTrials
        RandomSeed
    end
    methods
        function obj = DynamicForaging(varargin) % constructor            
            p = inputParser;
            
            % default parameters if none given
            p.addParameter('RewardProbabilities', [40 10]);
            p.addParameter('BlockLength', [50 100]);
            p.addParameter('MaxTrials', 1000);
            p.addParameter('RandomSeed', 1);
            p.parse(varargin{:});
            
            obj.RandomSeed = p.Results.RandomSeed;
            rng(obj.RandomSeed);
            obj.RewardProbabilities = p.Results.RewardProbabilities;
            obj.BlockLength = p.Results.BlockLength;
            obj.MaxTrials = p.Results.MaxTrials;
            obj.AllRewards = NaN(obj.MaxTrials, 2);
            obj.AllChoices = NaN(obj.MaxTrials, 2); % indexed as left and right
            obj.Trial = 0;
            obj.BlockSwitch = 1;            
            obj.BlockEnd = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockSwitch = [obj.BlockSwitch obj.BlockEnd];
        end
        function obj = inputChoice(obj, currChoice)
            % increment trial
            obj.Trial = obj.Trial + 1;
            
            % currChoice needs to be a 1x2 binary vector for [L R] choice
            if obj.NewBlock_Flag == true
                obj.BlockEnd = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
                obj.BlockSwitch = [obj.BlockSwitch obj.BlockEnd];
                obj.RewardProbabilities = obj.RewardProbabilities(2:-1:1); % switch reward probabilities
                obj.NewBlock_Flag = false;
            end
            
            % input choice
            if all(currChoice == [1 0]) % left choice
                obj.AllChoices(obj.Trial, :) = [1 0];
                if obj.BaitLeft_Flag == true % harvest L reward
                    obj.BaitLeft_Flag = false;
                    obj.AllRewards(obj.Trial, :) = [1 0];
                else
                    obj.AllRewards(obj.Trial, :) = [0 0];
                end
            elseif all(currChoice == [0 1]) % right choice
                obj.AllChoices(obj.Trial, :) = [0 1];
                if obj.BaitRight_Flag == true % harvest R reward
                    obj.BaitRight_Flag = false;
                    obj.AllRewards(obj.Trial, :) = [0 1];
                else
                    obj.AllRewards(obj.Trial, :) = [0 0];
                end
            else
                error('Choice should be [1 0] (left choice) or [0 1] (right choice)')
            end
            
            % generate rewards
            Lrand = randi([0 99]);
            Rrand = randi([0 99]);
            if Lrand < obj.RewardProbabilities(1) && obj.BaitLeft_Flag == false
                obj.BaitLeft_Flag = true;
            end
            if Rrand < obj.RewardProbabilities(2) && obj.BaitRight_Flag == false
                obj.BaitRight_Flag = true;
            end

            % next block
            if obj.Trial == obj.BlockSwitch(end)
                obj.NewBlock_Flag = true;
            end
        end
    end
end