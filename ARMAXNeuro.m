classdef ARMAXNeuro < handle
    %   ARMAXNeuro is a class for estimating intrinsic, seasonal, and exogenous
    %   memories and exogenous selectivity of neurons.
    %   
    %   To construct an instance use "ARMAXNeuo(neuronData)"
    %
    %   Copyright 2019 Mehran M. Spitmaan (mehran.m.spitman@gmail.com).
    %             2022 Ethan Trepka (ethan.trepka@gmail.com).
    
    properties (SetAccess = private)        
                
        timeIntervalsAbsolute           % Absolute time values for each bin and trial
               
        firingRateMat                   % Total firing rate for trial and bins

        firingRateMatAll                % Total firing rate for all trials and bins

        boundaries                       % Boundries for trials and bins
        
        totalDataPoint                  % Total number of datapoint
        
        firingRateMatMean               % Total Firing rate for bins (neural profile)
        
        firingRateMatMeanArranged       % Total Firing rate for bins (neural profile) rearranged for all data points
                
        fittingModelData                % Final Data needed for fitting the model
        
        display                         % Display messages in output
        
        paramSet                        % Data structure for parameter set and their boundaries
        
        fittingSet                      % Data structure for fitting procedure
        
        fittingResults                  % Data structure for fitting results
                        
        zscoring = 0;                   % Indicates whether zscoring the data
        
        timeMask                        % Time mask signals
        
        %% Data Structure for the Input Data
        neuronData
        % Contain all the input information we want to create the model
        %
        %   neuronData (Struct)
        %       .spikeTime:         All the spike times for entire experiment
        %
        %       .signalTime:        Structure contains the absolute time of all
        %                           signals for each trial.
        %                           Name of each field represents name of each
        %                           signal.
        %
        %       .signalValue:       Structure contains the values of all the possible
        %                           signals for each trial. There might be signals
        %                           with signalTime but without signalValue
        %                           Name of each field represents name of each
        %                           signal.
        %
        %       .endOfTrialSignal:  Name of the signal that indicates the end of the
        %                           trial. There might be a time distance
        %                           before/after "endOfTrialSignal" to the real end
        %                           of the trial which can be set by
        %                           "endOfTrialDist". (String)
        %
        %       .endOfTrialDist:    Time distanse from "endOfTrialSignal" that
        %                           indicated the real end of the trial. (msec)
        %
        %       .maxTrialLen:       Maximum time for each trial. Useful for
        %                           arranging data into bins. (msec)
        %                           The real value for trial length would be the
        %                           caculated by
        %                               "(signalTime.(endOfTrialSignal)(currTrial) +
        %                               endOfTrialDist) - ...
        %                               (signalTime.(endOfTrialSignal)(currTrial-1) +
        %                               endOfTrialDist);"
        %
        %       .beginOfTrialSignal:    Name of the signal that indicates the beginning of the
        %                               trial. (String)
        %
        %       .binSize:           Size of each bin (msec)
        %
        %       .trialRange:        Range of trials; Total number of trials.
        
        
        %% Data Structure for Model       
        
        component
        % Components contains all 3 possible memory types and exogenous
        % selectivity, each have their on data structure. 
        %
        %   .memShort (Struct for short memory)
        %       .defined:               Whether memory is defined or not
        %                               (Boolean)
        %
        %       .order:                 Indicates the memory order.
        %
        %       .varIdx:                Index of variables in total variable
        %                               vector.
        %
        %        
        %
        %   .memSeason (Struct for seasonal memory)
        %       .defined:               Whether memory is defined or not
        %                               (Boolean)
        %
        %       .order:                 Indicates the memory order
        %
        %       .varIdx:                Index of variables in total vriable
        %                               vector
        %
        %
        %   .memExo (Struct(s) for Exogenous memory)
        %       .name:                  Name of the exogenous memory
        %
        %           .signalName:        A set contains name(s) of the
        %                               signals involve in this memory       
        %
        %           .interaction:       Whether memory is based on
        %                               interaction of signals
        %
        %           .signalTime:        Name of the 'signalTime'
        %                               corresponding for this exo memory                                      
        %
        %
        %           .memOrder:          Indicates the memory order
        %
        %
        %   .exoSignal (Struct(s) for Exogenous signals)
        %       .name:                  Name of the exogenous memory
        %
        %           .signalName:        A set contains name(s) of the
        %                               signals involve in this memory       
        %
        %           .interaction:       Whether memory is based on
        %                               interaction of signals        
        %
        %           .effectiveTimeWin:  Time window that this memory is
        %                               effective (refer to timeMask function)
        %
    end
    
    methods
        % Constructor
        function obj = ARMAXNeuro(neuronData, display)
            % ARMAXNeuroLinear constructs an instance of this class
            %   You can introduce "neuronData" based on structure below:           
            obj.neuronData = neuronData;
            obj.display = display;
            
            % Initializing components
            % Short Memory
            component.memShort.defined = false;
            component.memShort.order = 0;
            
            % Seasonal Memory
            component.memSeason.defined = false;
            component.memSeason.order = 0;
                        
            % Exo Memory 
            component.memExo = struct;
            
            % Exo Signal
            component.exoSignal = struct; 
            
            obj.component = component;
        end
        
        % Initializing Data
        function initData(obj)
            %   Initiates the data and make it ready for further
            %   computation: 
            %   - Computing trial_alignment based on signals
            %   - Computing the firing rates
            %   - Computing ITI
            
            obj.displayMsg('Start initiation process...');
            
            % Compute absolute time intervals
            startPointEst = -obj.neuronData.maxTimeBeforeAlign;
            endPointEst = obj.neuronData.maxTimeAfterAlign;
            
            alignSignTime = obj.neuronData.signalTime.(obj.neuronData.alignTrialSignal);
            num_trials = size(alignSignTime,1);
            
            bin_range = [startPointEst:obj.neuronData.binSize:endPointEst];
            num_bins = size(bin_range,2);

            obj.timeIntervalsAbsolute = repmat(bin_range,[num_trials,1]) + ...
                repmat(alignSignTime,[1, num_bins]);
            
            % Calculating Firing Rate Matrix
            obj.firingRateMat = [];
            obj.firingRateMatAll = [];
            
            trial_starts = obj.neuronData.signalTime.(obj.neuronData.beginOfTrialSignal);
            trial_ends = obj.neuronData.signalTime.(obj.neuronData.endOfTrialSignal);
                
            for cntTrial = 1:num_trials              
                % All spike in current trial
                allSpikes = obj.neuronData.spikeTime;
                
                trialtimeIntervalsAbsolute = obj.timeIntervalsAbsolute(cntTrial,:);
                
                spikeCount = histc(allSpikes,trialtimeIntervalsAbsolute);
                spikeCount = spikeCount(1:end-1);
                
                if isempty(spikeCount)
                    spikeCount = zeros(length(trialtimeIntervalsAbsolute)-1,1);
                end
                
                if size(spikeCount,1)>1
                    spikeCount = spikeCount';
                end
                
                % store firing rate, including previous trial values
                obj.firingRateMatAll = [obj.firingRateMatAll; spikeCount/(obj.neuronData.binSize/1000)];
                
                % make the spike count equal to nan, if the bin is
                % before the start of or after the end of the current trial
                trial_start = trial_starts(cntTrial);
                trial_end = trial_ends(cntTrial);
                if cntTrial == 1
                    spikeCount((trialtimeIntervalsAbsolute+obj.neuronData.binSize) < trial_start) = nan;
                else
                    spikeCount((trialtimeIntervalsAbsolute-obj.neuronData.binSize) < trial_ends(cntTrial-1)) = nan;
                end
                spikeCount((trialtimeIntervalsAbsolute-obj.neuronData.binSize) > trial_end) = nan;
                
                % storing firing rate
                obj.firingRateMat = [obj.firingRateMat; spikeCount/(obj.neuronData.binSize/1000)];                
            end                     
            
           
            % calculate mean firing pattern in profile
            obj.firingRateMatMean = nanmean(obj.firingRateMat,1);                        

            % calculating boundaries
            obj.calcTrialBinBoundaries();
            
            % calculating model components
            obj.initMemShort();
            obj.initMemSeason();
            obj.initMemExo();
            obj.initExoSignal();

            obj.displayMsg('Initiation process is done!');
        end
               
        % Setup fitting procedure information
        function setFittingSet(obj, crossPer, crossIter, specificModel)     
            % crossPer - Cross-validation test data percentage
            % crossIter - Cross-validation instances 

            obj.displayMsg('Setup fitting parameters...');
            
            obj.fittingSet.crossPer = crossPer;
            obj.fittingSet.crossIter = crossIter;            
            obj.fittingSet.specificModel = specificModel;    
            
            obj.fittingSet.selectedComponents.seqOrderName = {}; % a dictonary for converting seq order num to name of components
            obj.fittingSet.selectedComponents.totalPermutationNumber = []; % total number of model permutations
            obj.fittingSet.selectedComponents.currentSelectedSeq = []; % The current selected sequence of components (binary seq)
            
            % Selection of ExoSignals
            obj.fittingSet.selectedComponents.exoSignal.selected = ~isempty(fieldnames(obj.component.exoSignal));
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'exoSignal';
            obj.fittingSet.selectedComponents.exoSignal.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            
            % Selection of short and seasonal memories
            obj.fittingSet.selectedComponents.memShort.selected = obj.component.memShort.defined;
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'memShort';
            obj.fittingSet.selectedComponents.memShort.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            
            obj.fittingSet.selectedComponents.memSeason.selected = obj.component.memSeason.defined;
            obj.fittingSet.selectedComponents.seqOrderName{end+1} = 'memSeason';
            obj.fittingSet.selectedComponents.memSeason.seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
            
            % Calculating exogenous memory selection
            if ~isempty(fieldnames(obj.component.memExo))
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                
                
                for cntExoMem = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntExoMem};
                    obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).selected = 1;
                    obj.fittingSet.selectedComponents.seqOrderName{end+1} = ['memExo.',exoMemNameTemp];
                    obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber = length(obj.fittingSet.selectedComponents.seqOrderName);
                end
            end
                        
            obj.fittingSet.selectedComponents.totalPermutationNumber = 2.^length(obj.fittingSet.selectedComponents.seqOrderName);
            obj.fittingSet.selectedComponents.currentSelectedSeq = ones(1,length(obj.fittingSet.selectedComponents.seqOrderName));
        end
        
        % Fit the model
        function fittingResults = fit(obj)                   
            obj.displayMsg('Running Specific Model...');
            cntModel = obj.fittingSet.specificModel-1;
            obj.fittingSet.currentModel = cntModel+1;

            % Set current selection Mask
            obj.fittingSet.selectedComponents.currentSelectedSeq = dec2bin(cntModel,length(obj.fittingSet.selectedComponents.currentSelectedSeq));
            obj.fittingSet.selectedComponents.currentSelectedSeq = obj.fittingSet.selectedComponents.currentSelectedSeq(end:-1:1);

            % Setup parameter set and boundaries
            obj.calcParamSet();

            % Fit a single model
            obj.fit_singleModel();
            fittingResults = obj.fittingResults.models;
            obj.displayMsg('End of Specific Model...');
        end
        
        % Fit a single model
        function fit_singleModel(obj)
            [X, y, X_time] = obj.getXandY();
            [vars, stat, ret] = obj.fittingFunction(X,y, X_time, obj.fittingSet.crossPer, obj.fittingSet.crossIter, obj.fittingSet.idxTrain, obj.fittingSet.idxTest, 0);
            
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.vars = vars;
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.stat = stat;
            obj.fittingResults.models{obj.fittingSet.currentModel}.res.ret = ret;
            
            obj.fittingResults.models{obj.fittingSet.currentModel}.paramSet = obj.paramSet;
        end
    end
    
    methods(Access = private)
        % Display a message
        function displayMsg(obj, msg)
            if obj.display
                disp(msg);
            end
        end     
        
        % Compute short term memory
        function initMemShort(obj)
            if (obj.component.memShort.defined)
                obj.displayMsg('Setup the neural data for Short memory...');
                
                dataPointIdx = 0;
                obj.fittingModelData.memShort.Rate = zeros(obj.totalDataPoint, obj.component.memShort.order);
                
                for cntTrials = obj.boundaries.trial
                    for cntBins = obj.boundaries.bin
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))            
                            % calculating the sliding window
                            binsSlidingBox = [cntBins-obj.component.memShort.order:cntBins-1];
                            
                            % Initializing data storages
                            ARSTemp = zeros(1,obj.component.memShort.order);
                            
                            cntId = 1;
                            for cntSlidingBox = binsSlidingBox
                                if cntSlidingBox > 0
                                    ARSTemp(cntId) = obj.firingRateMatAll(cntTrials,cntSlidingBox);
                                else
                                    ARSTemp(cntId) = 0;
                                end
                                cntId = cntId + 1;
                            end
                            
                            % Inverse the vectors
                            ARSTemp = ARSTemp(end:-1:1);
                            
                            dataPointIdx = dataPointIdx + 1;
                            obj.fittingModelData.memShort.Rate(dataPointIdx,:) = ARSTemp;                            
                        end
                    end
                end
            end
        end       
        
        % Compute seasonal memory
        function initMemSeason(obj)
            if (obj.component.memSeason.defined)
                obj.displayMsg('Setup the neural data for Seasonal memory...');
                
                dataPointIdx = 0;
                obj.fittingModelData.memSeason.Rate = zeros(obj.totalDataPoint, obj.component.memSeason.order);
                obj.fittingModelData.memSeason.Time = zeros(obj.totalDataPoint, obj.component.memSeason.order);

                for cntTrials = obj.boundaries.trial
                    for cntBins = obj.boundaries.bin
                        if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                            ARLTemp = obj.firingRateMatAll(cntTrials-obj.component.memSeason.order:cntTrials-1,cntBins);
                            timeL = obj.timeIntervalsAbsolute(cntTrials,cntBins) -...
                                obj.timeIntervalsAbsolute(cntTrials-obj.component.memSeason.order:cntTrials-1,cntBins);

                            % Inverse the vectors
                            ARLTemp = ARLTemp(end:-1:1);
                            timeL = timeL(end:-1:1);

                            dataPointIdx = dataPointIdx + 1;
                            obj.fittingModelData.memSeason.Rate(dataPointIdx,:) = ARLTemp;
                            obj.fittingModelData.memSeason.Time(dataPointIdx,:) = timeL;
                        end
                    end
                end
            end
        end
        
        % Compute exogenous memories
        function initMemExo(obj)
            if ~isempty(fieldnames(obj.component.memExo))
                obj.displayMsg('Setup the neural data for exo memories...');
                
                exoMemNames = fieldnames(obj.component.memExo);
                
                for cntExoMem = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntExoMem};
                    exoMemTemp = obj.component.memExo.(exoMemNameTemp);
                    
                    dataPointIdx = 0;
                    obj.fittingModelData.memExo.(exoMemNameTemp).Signal = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                    obj.fittingModelData.memExo.(exoMemNameTemp).Time = zeros(obj.totalDataPoint, exoMemTemp.memOrder);
                    
                    for cntTrials = obj.boundaries.trial
                        for cntBins = obj.boundaries.bin
                            if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                if ~exoMemTemp.interaction
                                    exoSigTemp = obj.neuronData.signalValue.(exoMemTemp.signalName)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                    exoTimeTemp = obj.timeIntervalsAbsolute(cntTrials,cntBins) -...
                                        obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                else
                                    for cntSignals = 1:length(exoMemTemp.signalName)
                                        exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoMemTemp.signalName{cntSignals})(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                    end
                                    
                                    exoSigTemp = obj.interactionCalc(exoSigTemp);
                                    exoTimeTemp = obj.timeIntervalsAbsolute(cntTrials,cntBins) -...
                                        obj.neuronData.signalTime.(exoMemTemp.signalTime)(cntTrials-exoMemTemp.memOrder:cntTrials-1);
                                end
                                
                                % Inverse the vectors
                                exoSigTemp = exoSigTemp(end:-1:1);
                                exoTimeTemp = exoTimeTemp(end:-1:1);
                                
                                dataPointIdx = dataPointIdx + 1;
                                obj.fittingModelData.memExo.(exoMemNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                                obj.fittingModelData.memExo.(exoMemNameTemp).Time(dataPointIdx,:) = exoTimeTemp;
                            end
                        end
                    end                    
                end
            end
        end
        
        % Compute exogenous signals
        function initExoSignal(obj)
            if ~isempty(fieldnames(obj.component.exoSignal))
                obj.displayMsg('Setup the neural data for exo signals...');
                
                exoSignalNames = fieldnames(obj.component.exoSignal);
                
                for cntExoSignal = 1:length(exoSignalNames)
                    exoSignalNameTemp = exoSignalNames{cntExoSignal};
                    exoSignalTemp = obj.component.exoSignal.(exoSignalNameTemp);
                    
                    dataPointIdx = 0;
                    obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal = zeros(obj.totalDataPoint, 1);
                    
                    for cntTrials = obj.boundaries.trial
                        for cntBins = obj.boundaries.bin
                            if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                                if ~exoSignalTemp.interaction
                                    % current time bigger than the signal
                                    if obj.timeIntervalsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoSignalTemp.signalName)(cntTrials))
                                        exoSigTemp = obj.neuronData.signalValue.(exoSignalTemp.signalName)(cntTrials);
                                    else % current time smaller than the signal
                                        exoSigTemp = 0; % mask anything smaller than signal
                                    end
                                else
                                    for cntSignals = 1:length(exoSignalTemp.signalName)
                                        % current time bigger than the signal
                                        if obj.timeIntervalsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(exoSignalTemp.signalName{cntSignals})(cntTrials))
                                            exoSigTemp(:,cntSignals) = obj.neuronData.signalValue.(exoSignalTemp.signalName{cntSignals})(cntTrials);
                                        else % current time smaller than the signal
                                            exoSigTemp(:,cntSignals) = 0; % mask anything smaller than signal
                                        end
                                    end
                                    
                                    exoSigTemp = obj.interactionCalc(exoSigTemp);
                                end
                                
                                % Inverse the vectors
                                exoSigTemp = exoSigTemp(end:-1:1);
                                dataPointIdx = dataPointIdx + 1;
                                obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal(dataPointIdx,:) = exoSigTemp;
                            end
                        end
                    end
                    
                    obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask = ones(obj.totalDataPoint, 1);
                    if ~isempty(obj.component.exoSignal.(exoSignalNameTemp).effectiveTimeWin)
                        obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask = obj.clcTimWin(obj.component.exoSignal.(exoSignalNameTemp).effectiveTimeWin, obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask);
                    end
                end
            end
        end       
        
        % Compute trial boundaries
        function calcTrialBinBoundaries(obj)
            % calculate maximum order
            allOrders(1) = obj.component.memShort.order;
            allOrders(2) = obj.component.memSeason.order;
            
            exoMemNames = fieldnames(obj.component.memExo);
            for cntExoMems = 1:length(exoMemNames)
                allOrders(end+1) = obj.component.memExo.(exoMemNames{cntExoMems}).memOrder;
            end
            
            maxOrder = max(allOrders);
            
            % calculating boundaries
            obj.boundaries.trial = maxOrder+1:size(obj.firingRateMat,1);
            obj.boundaries.bin = 1:size(obj.firingRateMat,2);
            
            % Calculating total datapoints
            obj.totalDataPoint = 0;
            obj.fittingModelData.output = [];
            obj.firingRateMatMeanArranged = [];
            for cntTrials = obj.boundaries.trial
                for cntBins = obj.boundaries.bin
                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                        obj.totalDataPoint = obj.totalDataPoint + 1;
                        obj.fittingModelData.output(end+1) = obj.firingRateMat(cntTrials,cntBins);
                        obj.firingRateMatMeanArranged(end+1) = obj.firingRateMatMean(cntBins);
                    end
                end
            end
        end
        
        % Compute time windows for time mask
        function tempMask = clcTimWin(obj, timeWin, tempMask)
            % initiate a temp mask
            dataPointIdx = 0;
    
            for cntTrials = obj.boundaries.trial
                for cntBins = obj.boundaries.bin
                    % clc
                    if ~isnan(obj.firingRateMat(cntTrials,cntBins))
                        dataPointIdx = dataPointIdx + 1;
                        % clc
                        % current time smaller than the mask begin time
                        if obj.timeIntervalsAbsolute(cntTrials,cntBins) < (obj.neuronData.signalTime.(timeWin.beginSig_Name)(cntTrials))
                            tempMask(dataPointIdx,:) = tempMask(dataPointIdx,:) * 0;
                        % current time biger than the mask end time
                        elseif obj.timeIntervalsAbsolute(cntTrials,cntBins) > (obj.neuronData.signalTime.(timeWin.endSig_Name)(cntTrials)) 
                            tempMask(dataPointIdx,:) = tempMask(dataPointIdx,:) * 0;
                        end

                    end
                end
            end
        end
        
        % Compute interaction data
        function data = interactionCalc(~, data)
            % Assumption is that size(data,2) shows the number of
            % parameters involved in interaction (presumably 2) and 
            % values can be either 1 or -1.
            
            data = prod(data,2);      
        end
        
        % Setup parameter set and boundaries
        function calcParamSet(obj)
            
            obj.displayMsg('Setup parameter settings...');
                        
            obj.paramSet.size.bias          = 1;
            obj.paramSet.size.exoSignal     = 1;
            obj.paramSet.size.memShort      = obj.component.memShort.order;
            obj.paramSet.size.memSeason     = obj.component.memSeason.order;
            
            if ~isempty(fieldnames(obj.component.memExo))
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                exoMemNameTemp = exoMemNames{1};
                exoSize = obj.component.memExo.(exoMemNameTemp).memOrder;
            else
                exoSize = 5;
            end
            
            obj.paramSet.size.memExo        = exoSize;
            
            obj.paramSet.flag.bias          = [];
            obj.paramSet.flag.exoSignal     = [];
            obj.paramSet.flag.memShort      = [];
            obj.paramSet.flag.memSeason     = [];
            obj.paramSet.flag.memExo        = [];
            
            obj.paramSet.effect.bias        = 0;
            obj.paramSet.effect.exoSignal   = 0;
            obj.paramSet.effect.memShort    = 0;
            obj.paramSet.effect.memSeason   = 0;
            obj.paramSet.effect.memExo      = 0;
            
            % Add bias term
            obj.paramSet.flag.bias = 1;
            currFlag = 1;
            obj.paramSet.effect.bias = 1;
            
            if ~isempty(fieldnames(obj.component.exoSignal)) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.exoSignal.seqNumber))
                obj.paramSet.flag.exoSignal = [currFlag+1:currFlag+obj.paramSet.size.exoSignal*(length(fieldnames(obj.component.exoSignal)))];
                currFlag = obj.paramSet.flag.exoSignal(end);
                obj.paramSet.effect.exoSignal = 1;
            end                        
            
            if (obj.component.memShort.defined) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memShort.seqNumber))
                obj.paramSet.flag.memShort = [currFlag+1:currFlag+obj.paramSet.size.memShort];
                currFlag = obj.paramSet.flag.memShort(end);
                obj.paramSet.effect.memShort = 1;
            end
            
            if (obj.component.memSeason.defined) && ...
                    str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memSeason.seqNumber))
                obj.paramSet.flag.memSeason = [currFlag+1:currFlag+obj.paramSet.size.memSeason];
                currFlag = obj.paramSet.flag.memSeason(end);
                obj.paramSet.effect.memSeason = 1;
            end
            
            if ~isempty(fieldnames(obj.component.memExo))
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                exoEffeciveNumTemp = 0; % store the numbers of effective exo memory components
                for cntN = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntN};
                    if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                        exoEffeciveNumTemp = exoEffeciveNumTemp + 1;
                    end
                end
                
                if exoEffeciveNumTemp>0
                    obj.paramSet.flag.memExo = [currFlag+1:currFlag+obj.paramSet.size.memExo*exoEffeciveNumTemp];
                    currFlag = obj.paramSet.flag.memExo(end);
                    obj.paramSet.effect.memExo = 1;
                end
            end
        end                
        
        % Compute output of the model
        function [X, y, X_time] = getXandY(obj)
            % predict firing rate in each bin minus mean firing rate in each
            % bin across trials
            
            y = obj.fittingModelData.output'-(obj.firingRateMatMeanArranged');
            
            % zscoring if needed
            if obj.zscoring == 1
                obj.zscoreData();
            else
                obj.fittingModelData.snapShot = obj.fittingModelData;
            end
            
            X = [];
            X_time = [];
            
            % calculating exogenous signal(s) component
            if obj.paramSet.effect.exoSignal
                exoSignalNames = fieldnames(obj.fittingModelData.exoSignal);
                for cntExoSignal = 1:length(exoSignalNames)
                    exoSignalNameTemp = exoSignalNames{cntExoSignal};
                    exoSignalPart(:,cntExoSignal)   = obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).mask .* obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).Signal;
                end
                exoSignalPart(isnan(exoSignalPart))     =0;
                X = [X, exoSignalPart];
            end
            
            % Calculating short memory component
            if obj.paramSet.effect.memShort
                memShortPart = obj.fittingModelData.snapShot.memShort.Rate;
                memShortPart(isnan(memShortPart))       = 0;
                X = [X, memShortPart];
            end
            
            % Calculating long memory component
            if obj.paramSet.effect.memSeason
                memSeasonPart     =  obj.fittingModelData.snapShot.memSeason.Rate;
                memSeasonPart(isnan(memSeasonPart))     = 0;
                X = [X, memSeasonPart];
            end
            
            % Calculating exogenous memory component
            if obj.paramSet.effect.memExo
                exoMemNames = fieldnames(obj.fittingModelData.memExo);
                cntExoMemTemp = 1;
                for cntExoMem = 1:length(exoMemNames)
                    exoMemNameTemp = exoMemNames{cntExoMem};
                    if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                        exoMemTemp = obj.fittingModelData.snapShot.memExo.(exoMemNameTemp);
                        memExoPart(:,cntExoMemTemp:cntExoMemTemp+size(exoMemTemp.Signal,2)-1) = exoMemTemp.Signal.*obj.firingRateMatMeanArranged';%.*exoMemTemp.Time/1000;
                        memExoTime(:,cntExoMemTemp:cntExoMemTemp+size(exoMemTemp.Signal,2)-1) = exoMemTemp.Time/1000;
                        %memExoPart(:,cntExoMemTemp+size(exoMemTemp.Signal,2):cntExoMemTemp+2*size(exoMemTemp.Signal,2)-1) = exoMemTemp.Signal.*obj.firingRateMatMeanArranged';
                        cntExoMemTemp = cntExoMemTemp + size(exoMemTemp.Signal,2);
                    end
                end
                memExoPart(isnan(memExoPart))           = 0;
                X = [X, memExoPart];
                X_time = [X_time, memExoTime];
            end            
        end
        
        % ZScoreData
        function zscoreData(obj)
            % Exo Signals
            exoSignalNames = fieldnames(obj.fittingModelData.exoSignal);
            
            for cntExoSignal = 1:length(exoSignalNames)
                exoSignalNameTemp = exoSignalNames{cntExoSignal};
                obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).Signal = nanzscore(obj.fittingModelData.exoSignal.(exoSignalNameTemp).Signal);     
                obj.fittingModelData.snapShot.exoSignal.(exoSignalNameTemp).mask = obj.fittingModelData.exoSignal.(exoSignalNameTemp).mask;
            end
            
            % Short Memory
            obj.fittingModelData.snapShot.memShort.Rate = nanzscore(obj.fittingModelData.memShort.Rate);
            
            % Seasonal Memory
            obj.fittingModelData.snapShot.memSeason.Rate = nanzscore(obj.fittingModelData.memSeason.Rate);
            
            % Exogenous Memory
            exoMemNames = fieldnames(obj.fittingModelData.memExo);
            
            for cntExoMem = 1:length(exoMemNames)
                exoMemNameTemp = exoMemNames{cntExoMem};
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).Signal = nanzscore(obj.fittingModelData.memExo.(exoMemNameTemp).Signal);
                obj.fittingModelData.snapShot.memExo.(exoMemNameTemp).Time = (obj.fittingModelData.memExo.(exoMemNameTemp).Time);

            end
        end

        function [vars, stat, ret] = fittingFunction(obj, X,y,X_time, crossPer, crossIter, idxTrain, idxTest, display)
        % A general fitting function for cross-validated linear regression.
        % 
        %   [USAGE]
        %   [vars, stat, ret] = fit_lienar_cv(X,y, crossPer, crossIter, display)
        %
        %   [INPUT]
        %   X:          Input predictor matrix.
        %   y:          output to predict
        %
        %   crossPer:   Percentage of test data considred for each cross-validation
        %               instance. [0:1]
        % 
        %   crossIter:  Number of cross-validation instances.
        %
        %   display:    Display iterations for cross-validation process 
        %               (true, false)
        %               
        %
        %
        %   [OUTPUT]
        %   vars:       Final solution.
        %   
        %   stat:       Containing information about the solution. 
        %               stat.all:     Stats of all data
        %               stat.train:     Stats of train data
        %               stat.test:      Stats of test data
        %                       .N:     Number of data (all, test, train)
        %                       .K:     Number of variables (all, test, train)
        %                       .RSS:   Residual sum of squares (all, test, train)
        %                       .AIC:   AIC of the data (all, test, train)

        %                       .P_test:P-values of of all parameters
        %
        %
        %   ret:        Extra return values
        %                       .fvalFinal      Best f-value
        %                       .fvarsFinal     Best solution
        %
        %                       .fvalTest       F-values for test data (all instances and iterations)
        %                       .fvalTrain      F-values for train data (all instances and iterations)
        %                       .fvarsAll       Vars (all instances and iterations)
        % 
        %                       .idxSortFinal   Best models in each instance
        % 
        %                       .fvalFinalset   F-values for test data (all instances)
        %                       .fvarsFinalset  Vars (all instances)
        %
        % Copyright 2019 Mehran M. Spitmaan (mehran.m.spitman@gmail.com); 
        %           2022 Ethan Trepka


        %% if bias only model, create mock X variable and change model spec
        bias_only = 0;
        if numel(X) == 0
            X = zeros(size(y));
            bias_only = 1;
        end

        %% Iteration over cross-validation instances
        mdlsTrain = {};
        for cntInstance = 1:crossIter
            X_train = X(idxTrain{cntInstance},:);
            y_train = y(idxTrain{cntInstance});
            
            X_test = X(idxTest{cntInstance},:);

            y_test = y(idxTest{cntInstance});
            
            if obj.paramSet.effect.memExo
                X_time_train = X_time(idxTrain{cntInstance},:);
                X_time_test = X_time(idxTest{cntInstance},:);
            end
            
            % Training Part
            if bias_only
                mdl = fitlm(X_train, y_train, 'constant');
            else
                if ~obj.paramSet.effect.memExo
                    mdl = fitlm(X_train, y_train);
                else          
                    % get names of exogenous memory components
                    exoMemNames = fieldnames(obj.fittingModelData.memExo);
                    num_mem = 0; % store the numbers of effective exo memory components
                    for cntN = 1:length(exoMemNames)
                        exoMemNameTemp = exoMemNames{cntN};
                        if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                            num_mem = num_mem + 1;
                        end
                    end
            
                    % get saved hyperparameters
                    tau_upper = obj.component.memExo.(exoMemNameTemp).tauLim(2);
                    tau_lower = obj.component.memExo.(exoMemNameTemp).tauLim(1);
                    amp_upper = obj.component.memExo.(exoMemNameTemp).ampLim(2);
                    amp_lower = obj.component.memExo.(exoMemNameTemp).ampLim(1);
                    logn_params = obj.component.memExo.(exoMemNameTemp).logn_params;
                    mem_order = obj.component.memExo.(exoMemNameTemp).memOrder;

                    % make 10 initial guesses at tau parameter, then fit
                    num_iter = 10;
                    
                    all_taus = zeros(num_iter, num_mem);
                    all_r2 = zeros(num_iter,1);

                    X_fit = [ones(size(X_train,1),1), X_train(:,1:end-num_mem*mem_order), zeros(size(X_train,1),num_mem)];
                    
                    best_r2 = 0;
                    for i = 1:num_iter
                        taus = lognrnd(logn_params(1), logn_params(2), 1,num_mem);
                        while sum(taus>tau_upper)>0 | sum(taus<tau_lower)>0
                           taus = lognrnd(logn_params(1), logn_params(2), 1,num_mem);
                        end
                      
                        X_fit = obj.calcExoDesign(X_train, X_time_train, taus);
                        mdl_new = fitlm(X_fit(:,2:end), y_train);

                        if mdl_new.Rsquared.Ordinary>best_r2 || best_r2 == 0
                            best_r2 = mdl_new.Rsquared.Ordinary;
                            mdl_best = mdl_new;
                        end
                    end
                    init_params = [mdl_best.Coefficients.Estimate', taus];
                    
                    options             = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective','MaxIter',5000, 'MaxFunEvals',10000, 'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1e-4);%,'CheckGradients', true);
                    % options             = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','MaxIter',5000, 'MaxFunEvals',10000);
                    
                    options.Display     = 'off';
                    minFuncTrain                        = @(v) obj.minFunc(v, num_mem, mem_order, y_train, X_train, X_time_train);
                    lowerLim = [repmat(-inf, 1, size(X_fit,2)-num_mem),repmat(amp_lower, 1, num_mem),repmat(tau_lower, 1, num_mem)];
                    upperLim = [repmat(inf, 1, size(X_fit,2)-num_mem),repmat(amp_upper, 1, num_mem),repmat(tau_upper, 1, num_mem)];
                    [solution,RESNORML,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] =...
                        lsqnonlin(minFuncTrain,init_params,lowerLim,upperLim,options);
                    taus = solution(end-num_mem+1:end)';
                    
                    X_fit = obj.calcExoDesign(X_train, X_time_train, taus);
                    
                    mdl_new = fitlm(X_fit(:,2:end), y_train);
                    
                    X_train = obj.calcExoDesign(X_train, X_time_train, taus);
                    X_test = obj.calcExoDesign(X_test, X_time_test, taus);
                    
                    mdl = fitlm(X_train, y_train);
                end
            end

            mdlsTrain{cntInstance} = mdl;
            varsTrain(:,cntInstance)    = mdl.Coefficients.Estimate;
            
            if obj.paramSet.effect.memExo
                tausTrain(:,cntInstance) = taus;
            end
            fvalTrain(cntInstance)    = mdl.SSR;

            test_pred = mdl.predict(X_test);
            test_error = sum((y_test-test_pred).^2);

            fvalTest(cntInstance)     = test_error;

            if display
                [cntInstance]
            end    
        end
        
        if obj.paramSet.effect.memExo
            varsTrain(end+1:end+num_mem,:) = tausTrain;
        end
        %% calculating best model based on test SSE
        [~, idxTemp]        = min(fvalTest);

        fvalTestBest        = fvalTest(idxTemp);
        varsTrainBest       = varsTrain(:,idxTemp);
        mdlsTrainBest       = mdlsTrain{idxTemp};

        % Calculating median of best values over cross-validation instances
        fvalFinal           = median(fvalTest);
        fvarsFinal          = varsTrainBest;

        % setting return values
        ret.fvalFinal       = fvalFinal;        % best f-value
        ret.fvarsFinal      = fvarsFinal;       % best solution

        ret.fvalTest        = fvalTest;         % f-values for test data (all instances and iterations)
        ret.fvalTrain       = fvalTrain;        % f-values for train data (all instances and iterations)
        ret.fvarsAll        = varsTrain;        % vars (all instances and iterations)

        ret.idxSortFinal    = idxTemp;          % best models in each instance

        ret.fvalFinalset    = fvalTestBest;     % f-values for test data (all instances)
        ret.fvarsFinalset   = varsTrainBest;    % vars (all instances)

        % setting final return value
        vars                = fvarsFinal';

        %% Calculating STATS for TRAIN dataset
        statTemp = [];
        for cntInstance = 1:crossIter

            X_train = X(idxTrain{cntInstance},:);
            y_train = y(idxTrain{cntInstance});

            mdl                            = mdlsTrain{cntInstance};
            statTemp{cntInstance}.RSS           = mdl.SSE;   
            statTemp{cntInstance}.N             = length(y_train);
            statTemp{cntInstance}.K             = size(varsTrain,1);

            statTemp{cntInstance}.AIC           = 2*statTemp{cntInstance}.K + statTemp{cntInstance}.N * log(statTemp{cntInstance}.RSS);               

        end

        stat.train = statTemp;


        %% Calculating STATS for TEST dataset
        statTemp = [];
        for cntInstance = 1:crossIter
            X_test = X(idxTest{cntInstance},:);
            if obj.paramSet.effect.memExo
                X_test_time = X_time(idxTest{cntInstance},:);
                taus = varsTrain(end-num_mem+1:end,cntInstance);
                X_test = obj.calcExoDesign(X_test, X_time_test, taus);
            end      
            y_test = y(idxTest{cntInstance});
            
            mdl                            = mdlsTrain{cntInstance};
            test_pred = mdl.predict(X_test);
            test_error = sum((y_test-test_pred).^2);

            statTemp{cntInstance}.RSS           = test_error;   
            statTemp{cntInstance}.N             = length(y_test);
            statTemp{cntInstance}.K             = size(varsTrain,1);

            statTemp{cntInstance}.AIC           = 2*statTemp{cntInstance}.K + statTemp{cntInstance}.N * log(statTemp{cntInstance}.RSS);               
        end

        stat.test = statTemp;



        %% Calculating STATS for ALL dataset
        if obj.paramSet.effect.memExo
           taus = varsTrainBest(end-num_mem+1:end);
           X = obj.calcExoDesign(X, X_time, taus);
        end      
        
        statTemp = [];

        mdl                            = mdlsTrainBest;
        test_pred = mdl.predict(X);
        test_error = sum((y-test_pred).^2);

        test_sst = sum((y-mean(y)).^2);

        statTemp.RSS           = test_error;   
        statTemp.Rsquare       = 1-(test_error/test_sst);
        statTemp.N             = length(y);
        statTemp.K             = size(varsTrain,1);
        statTemp.AIC           = 2*statTemp.K + statTemp.N * log(statTemp.RSS);               

        statTemp.P_test        = mdlsTrainBest.Coefficients.pValue;

        stat.all = statTemp;

        end

        function X = calcExoDesign(obj, X, X_time, taus)
            % get names of exogenous memory components
            exoMemNames = fieldnames(obj.fittingModelData.memExo);
            num_mem = 0; % store the numbers of effective exo memory components
            for cntN = 1:length(exoMemNames)
                exoMemNameTemp = exoMemNames{cntN};
                if str2num(obj.fittingSet.selectedComponents.currentSelectedSeq(obj.fittingSet.selectedComponents.memExo.(exoMemNameTemp).seqNumber))
                    num_mem = num_mem + 1;
                end
            end
            mem_size = round(size(X_time,2)/num_mem);

            % get part of design matrix that is unmodified
            X_base = X(:,1:end-size(X_time,2));
            X_exo = zeros(size(X,1),num_mem);
            
            % ge part of matrix that will be modified
            X_init = X(:,end-size(X_time,2)+1:end);
            
            % setup final matrix
            X = [X_base, X_exo];

            % update final matrix
            cntMemExo = 1;
            for j = 1:num_mem
                X_exo(:,j) = nansum(X_init(:,cntMemExo:cntMemExo+mem_size-1).*exp(-X_time(:,cntMemExo:cntMemExo+mem_size-1)/taus(j)),2);
                cntMemExo = cntMemExo + mem_size;
            end

            X(:,end-num_mem+1:end) = X_exo;
        end
        
        function [loss, J] = minFunc(obj, params, num_mem, mem_order, y_train, X_train, X_time_train)
            % get exogenous memory timescales
            taus = params(end-num_mem+1:end);
               
            % get linear parameters
            linear_pars = params(1:end-num_mem);
            
            % compute new design matrix based on taus
            X_fit = obj.calcExoDesign(X_train, X_time_train, taus);

            % generate model prediction
            y_pred = X_fit*linear_pars';
            
            % generate loss for lsqonlin
            loss = y_train-y_pred;
            
            % compute Jacobian
            if nargout > 1
                J = zeros(length(y_pred),length(params));
                % for all linear parameters the Jacobian = -X
                J(:,1:length(linear_pars)) = -X_fit(:,1:length(linear_pars));
                
                % for tau parameters, Jacobian is more complex and is
                % computed as shown below, correctness can be verified using 'grad
                % check' operation
                si = size(X_train,2)-size(X_time_train,2);
                cntMemExo = 1;
                mem_size = mem_order;
                for j = 1:num_mem
                    amp = linear_pars(end-num_mem+j);
                    pt1 = X_train(:,si+cntMemExo:si+cntMemExo+mem_size-1).*amp.*X_time_train(:,cntMemExo:cntMemExo+mem_size-1)/(taus(j)^2);
                    pt2 = exp(-X_time_train(:,cntMemExo:cntMemExo+mem_size-1)/taus(j));
                    J(:,length(linear_pars)+j) = -sum(pt1.*pt2,2);
                    cntMemExo = cntMemExo+mem_size;
                end
            end
        end
    end
    
    methods(Access = public)
        function timeMaskSignal(obj, name, beginSig_Name, endSig_Name)
            % creates timeMask (tm) object based on signals
            obj.timeMask.(name).beginSig_Name = beginSig_Name;
            obj.timeMask.(name).endSig_Name = endSig_Name;
        end
        
        function addMemShort(obj, order, interval)
            %   Initialize data structure for short memory
            %       .order:                 Indicates the memory order
            obj.component.memShort.defined = true;
            obj.component.memShort.order = order;
            obj.component.memShort.interval = interval;
        end
        
        function addMemSeason(obj, order, interval)
            %   Initialize data structure for seasonal memory
            %       .order:                 Indicates the memory order
            obj.component.memSeason.defined = true;
            obj.component.memSeason.order = order;
            obj.component.memSeason.interval = interval;
        end
                        
        function addMemExo(obj, name, signalNames, interaction, signalTime, memOrder, ampLim, tauLim, lognParams)
            %   Add and initialize data structure for exo memory
            %
            %           .name:                  Name of the exogenous memory
            %
            %           .signalName:        A set contains name(s) of the
            %                               signals involve in this memory       
            %
            %           .interaction:       Whether memory is based on
            %                               interaction of signals
            %
            %           .signalTime:        Name of the 'signalTime'
            %                               corresponding for this exo memory                                      
            %
            %           .memOrder:          Indicates the memory order 
            %           .ampLim:            Amplitude parameter limits
            %                               for fitting: amp_lower, amp_upper
            %           .tauLim:            Timescale parameter limits
            %                               for fitting: tau_lower, tau_upper
            
            obj.component.memExo.(name).signalName = signalNames;
            obj.component.memExo.(name).interaction = interaction;
            obj.component.memExo.(name).signalTime = signalTime;
            obj.component.memExo.(name).memOrder = memOrder;
            obj.component.memExo.(name).ampLim = ampLim;
            obj.component.memExo.(name).tauLim = tauLim;
            obj.component.memExo.(name).logn_params = lognParams;
        end
                        
        function addExoSignal(obj, name, signalNames, interaction, effectiveTimeWin)
            %   Add and initialize data structure for exo signal
            %       .name:                  Name of the exogenous memory
            %
            %           .signalName:        A set contains name(s) of the
            %                               signals involve in this memory
            %
            %           .interaction:       Whether memory is based on
            %                               interaction of signals
            %
            %           .effectiveTimeWin:  Time window that this memory is
            %                               effective (refer to timeMask function)

            obj.component.exoSignal.(name).signalName = signalNames;
            obj.component.exoSignal.(name).interaction = interaction;
            obj.component.exoSignal.(name).effectiveTimeWin = obj.timeMask.(effectiveTimeWin);
        end
    
        function [fitting_results, model_selection, output] = fitAllModelSelection(obj, num_models, crossIter)
            fitting_results = obj.fitAll(num_models, crossIter);
            model_selection = obj.doModelSelection(fitting_results, num_models);
            output = obj.constructOutput(fitting_results, model_selection);
        end
        
        function fitting_results = fitAll(obj, num_models, crossIter)
            % Setup cross-validation instances, use same instances for all models
            crossPer = 1/crossIter;
            
            totNum              = length(obj.firingRateMatMeanArranged); % Total number of samples
            foldNum             = floor(totNum * crossPer); % Number of samples in each fold

            idxTemp = [];
            for cntIter             = 1:crossIter
                rndtemp             = randperm(totNum);
                idxTemp{cntIter}    = rndtemp(1:foldNum);
            end

            idxTrain            = []; % Indices of train dataset
            idxTest             = []; % Indices of test dataset
            for cntIter = 1:crossIter    
                idxTrain{cntIter}   = setdiff([1:totNum],idxTemp{cntIter});
                idxTest{cntIter}    = [idxTemp{cntIter}];
            end

            obj.fittingSet.idxTrain = idxTrain;
            obj.fittingSet.idxTest = idxTest;
            
            % fit all models
            fitting_results = {};

            for specificModel = 1:num_models
                fit_stream = getByteStreamFromArray(obj);
                Model = getArrayFromByteStream(fit_stream);
                Model.setFittingSet(crossPer, crossIter, specificModel);
                model_results = Model.fit();
                fitting_results{specificModel} = model_results{specificModel};
            end
        end
        
        function model_selection = doModelSelection(obj, fitting_results, num_models)
            % variables to fill for a single neuron (for model selection)
            rssSetModel = nan(1,num_models);
            varsSetModel = {};
            pvalSetModel = {};
            flagSetModel = {};
            r2SetModel = nan(1,num_models);
            
            
            % iterate over all models
            rssSets = [];
            for cntModel = 1:num_models
                varsSetModel{cntModel} = fitting_results{1,cntModel}.res.vars;
                pvalSetModel{cntModel} = fitting_results{1,cntModel}.res.stat.all.P_test;
                r2SetModel(cntModel) = fitting_results{1,cntModel}.res.stat.all.Rsquare;
                rssSetTemp = [];
                
                for cntIns = 1:length(fitting_results{1,cntModel}.res.stat.test)
                    rssSetTemp(cntIns) = fitting_results{1,cntModel}.res.stat.test{1,cntIns}.RSS;
                end
                
                % rssSetModel is the average residual sums of sqaures for
                % cross-validation instances
                rssSetModel(cntModel) = mean(rssSetTemp);
                rssSets = [rssSets; rssSetTemp];
                % flagSetModel describes the components included in the model
                flagSetModel{cntModel} = fitting_results{1,cntModel}.paramSet.flag;
            end
            
            % perform  model selection, i.e., select the model with minimum
            % residual sums of squares (VT is RSS value, IT is model index)
            [VT, IT] = min(rssSetModel);
            
            % translate index into binary for components included in the model
            min_digits = log(num_models)/log(2);
            
            compartmentSeq = dec2bin(IT-1,min_digits);
            compartmentSeq = compartmentSeq(end:-1:1);
            compartmentSeq = compartmentSeq(4:end);
            
            % get parameters, pvalues, and flags for the best model for this neuron
            varsSetBest = varsSetModel{IT};
            pvalSetBest = pvalSetModel{IT};
            flagSetBest = flagSetModel{IT};
            flagSetAll  = flagSetModel{num_models};
            
            % stores the values for all models in allNeuronSt
            model_selection.rssSetModel = rssSetModel;
            model_selection.varsSetModel = varsSetModel;
            model_selection.pvalSetModel = pvalSetModel;
            model_selection.r2SetBest = r2SetModel(IT);
            model_selection.r2SetModel = r2SetModel;
            model_selection.IT = IT;
            model_selection.VT = VT;
            model_selection.flagSetModel = flagSetModel;
            
            % store the RSS for the best model
            model_selection.rssSet = VT;
            
            % Fields
            cellNames = fieldnames(flagSetBest);
            for cntCell = [1:4]
                cellNameTemp = cellNames{cntCell};
                if ~isempty(flagSetBest.(cellNameTemp))
                    memExoTempIdx_General = flagSetAll.(cellNameTemp);
                    memExoTempIdx_Special = flagSetBest.(cellNameTemp);
                    
                    if ~strcmp(cellNameTemp,'mean')
                       model_selection.varsSet(memExoTempIdx_General) =...
                            varsSetBest(memExoTempIdx_Special);
                        
                       model_selection.pvalSet(memExoTempIdx_General) =...
                            pvalSetBest(memExoTempIdx_Special);
                    end
                end
            end
            cellNameTemp = 'memExo';
            memExoSize = 2;
            
            % first, count number of exogenous components in best model
            exoEffeciveNumTemp = 0; % store the numbers of effective exo memory components
            
            for cntN = 1:length(compartmentSeq)
                exoMemTemp = compartmentSeq(cntN);
                if str2num(exoMemTemp)
                    exoEffeciveNumTemp = exoEffeciveNumTemp + 1;    
                end
            end
            num_mem = exoEffeciveNumTemp;
            num_mem_all = length(fieldnames(obj.component.memExo));
            % next, add components to set 
            exoEffectiveNumTemp = 1; % store the numbers of effective exo memory components

            for cntN = 1:length(compartmentSeq)
                % gets 1 or 0 from binary string created earlier which indicates
                % whether or not particular exogenous component was included
                exoMemTemp = compartmentSeq(cntN);
                memExoTempIdx_General = [cntN, cntN+num_mem_all];
                memExoTempIdx_General = flagSetAll.(cellNameTemp)(memExoTempIdx_General);
                if str2num(exoMemTemp)
                    memExoTempIdx_Special = [exoEffectiveNumTemp, exoEffectiveNumTemp+num_mem];

                    memExoTempIdx_Special = flagSetBest.(cellNameTemp)(memExoTempIdx_Special);

                    model_selection.varsSet(memExoTempIdx_General) =...
                        varsSetBest(memExoTempIdx_Special);

                    model_selection.pvalSet(memExoTempIdx_General) = pvalSetBest(memExoTempIdx_Special(1));
                    
                    exoEffectiveNumTemp = exoEffectiveNumTemp + 1;
                 else
                     model_selection.varsSet(memExoTempIdx_General) = 0;
                     model_selection.pvalSet(memExoTempIdx_General) = 0;
                end
            end
        end
        
        function output = constructOutput(obj, fitting_results, model_selection)
            % 1. get variable values, stored in ('var'), and p-values,
            % stored in ('p') for each parameter of interest
            % each parameter is assigned a p-value - for memory components components,
            % this is the minimum p-value of any of the memory terms
            flagSetAll = model_selection.flagSetModel{end};
            flagSetBest = model_selection.flagSetModel{model_selection.IT};


            biasIdx = flagSetAll.bias;
            exoSignalIdx = flagSetAll.exoSignal;
            shortIdx = flagSetAll.memShort;
            longIdx = flagSetAll.memSeason;
            memExoIdx = flagSetAll.memExo;
            
            varsSet = model_selection.varsSet;
            pvalSet = model_selection.pvalSet;
            
            % output struct
            output = struct;
            
            % bias term
            output.var.bias = varsSet(biasIdx);
            output.p.bias = pvalSet(biasIdx);
            
            % exoSignal
            exoSignalNames = fieldnames(obj.component.exoSignal);
            for cnt = 1:length(exoSignalNames)
                name = exoSignalNames{cnt};
                output.var.("exo_" + name) = varsSet(exoSignalIdx(cnt));
                output.p.("exo_" + name) = pvalSet(exoSignalIdx(cnt));
            end
            
            % memShort
            short_lag = obj.component.memShort.interval;
            short_coeffs = varsSet(shortIdx);
            short_tau = obj.getARTau(short_coeffs, short_lag);
            output.var.tau_intrinsic = short_tau;
            output.p.tau_intrinsic = min(pvalSet(shortIdx));

            % memSeason
            seasonal_lag = obj.component.memSeason.interval;
            seasonal_coeffs = varsSet(longIdx);
            seasonal_tau = obj.getARTau(seasonal_coeffs, seasonal_lag);
            output.var.tau_seasonal = seasonal_tau;
            output.p.tau_seasonal = min(pvalSet(longIdx));
            
            % memExo
            memExoNames = fieldnames(obj.component.memExo);
            num_mem = length(memExoNames);
            if length(memExoNames)>0
                sidx = memExoIdx(1);
            end
            for cnt = 1:length(memExoNames)
                name = memExoNames(cnt);
                output.var.("amp_" + name) = varsSet(sidx+cnt-1);
                output.var.("tau_" + name) = varsSet(sidx+cnt-1+num_mem);
                output.p.("tau_" + name) = pvalSet(sidx+cnt-1);
                output.p.("amp_" + name) = pvalSet(sidx+cnt-1);
            end
            
            fields = fieldnames(output.var);
            for i = 1:length(fields)
                f = fields{i};
                if output.var.(f) == 0
                    output.var.(f) = NaN;
                    output.p.(f) = NaN;
                end
            end
            
            % 2. get which components are included in all models, and the best fitting
            % model
            output.comp.labels = ["bias", "exogenous", "intrinsic", "seasonal"];
            for cnt = 1:length(memExoNames)
                output.comp.labels = [output.comp.labels, memExoNames(cnt)];
            end
            
            output.comp.included = [];
            allFields = fieldnames(flagSetAll);
            for cnt = 1:length(allFields)
                field = allFields{cnt};
                if ~strcmp(field, 'memExo')
                    if sum(flagSetBest.(field))>0
                        output.comp.included = [output.comp.included,1];
                    else
                        output.comp.included = [output.comp.included,0];
                    end
                else
                    for i = 1:length(memExoNames)
                        name = memExoNames(i);
                        if ~isnan(output.var.("tau_" + name))
                            output.comp.included = [output.comp.included,1];
                        else
                            output.comp.included = [output.comp.included,0];
                        end                          
                    end
                end
            end
            
            % 3. get R^2 value of best fitting model
            output.stat.r2_best = model_selection.r2SetBest;
            
            % 4. get R^2 difference between best model and model without
            % particular component           
            num_models = length(model_selection.r2SetModel);
            min_digits = log(num_models)/log(2);
            for model_num = 1:num_models
                compartmentSeq = dec2bin(model_num-1,min_digits);
                compartmentSeq = compartmentSeq(end:-1:1);
                %compartmentSeq = compartmentSeq(4:end);
                cs = str2num(compartmentSeq(:));
                if sum(cs) == length(cs)-1 % if this model has all but one component
                    delta_r2 = model_selection.r2SetModel(end)-model_selection.r2SetModel(model_num);
                    comp_idx = find(~cs);
                    label = output.comp.labels(comp_idx+1);
                    output.stat.delta_r2.(label) = delta_r2;
                end
            end
        end
        
        function tau = getARTau(obj, ar_coeffs, lag)
            n = length(ar_coeffs);
            eigMat = [ar_coeffs;[eye(n-1) zeros(n-1,1)]];
            lamdas = eig(eigMat);
            tau = abs(-lag./log(max(abs(lamdas))));
        end
    end
end

%% helper function for zscoring model data with nans
function [z,mu,sigma] = nanzscore(x,flag,dim)
%ZSCORE Standardized z score.
%   Z = ZSCORE(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-MEAN(X)) ./ STD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE(X) also returns MEAN(X) in MU and STD(X) in SIGMA.
%
%   [...] = ZSCORE(X,1) normalizes X using STD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  ZSCORE(X,0) is the same as
%   ZSCORE(X).
%
%   [...] = ZSCORE(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1
%   to use N.
%
%   See also MEAN, STD.

%   Copyright 1993-2015 The MathWorks, Inc. 


% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = x; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = nanmean(x,dim);
sigma = nanstd(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end
