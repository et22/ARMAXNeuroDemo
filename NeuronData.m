classdef NeuronData < handle
    %   NeuronData is a class for representing the data structure for
    %   trial-based neural data.
    %   
    %   To construct an instance use "NeuronData(spikeTime, binSize)". See
    %   testARMAXNeuro for a demonstration of its use.
    %
    %   Copyright 2022 Ethan Trepka (ethan.trepka@gmail.com).
    
   properties
      spikeTime
      
      binSize

      signalTime
      signalValue
      
      beginOfTrialSignal
      endOfTrialSignal
      
      alignTrialSignal
      maxTimeBeforeAlign
      maxTimeAfterAlign
   end
   methods
        % Constructor
        function obj = NeuronData(spikeTime, binSize)
            % NeuronData constructs an instance of this class
            obj.spikeTime = spikeTime;
            obj.binSize = binSize;
        end
   end
        
   methods(Access = public)        
        function addSignalTime(obj, name, time)           
            obj.signalTime.(name) = time;
        end
        
        function addSignalValue(obj, name, value)
            obj.signalValue.(name) = value;
        end
        
        function addTrialBeginEndSignal(obj, beginOfTrialSignal, endOfTrialSignal)
            obj.beginOfTrialSignal = beginOfTrialSignal;
            obj.endOfTrialSignal = endOfTrialSignal;
        end
        
        function addAlignTrialSignal(obj, signalName, maxTimeBeforeAlign, maxTimeAfterAlign)
            obj.alignTrialSignal = signalName;
            obj.maxTimeBeforeAlign = maxTimeBeforeAlign;
            obj.maxTimeAfterAlign = maxTimeAfterAlign;
        end
    end
end