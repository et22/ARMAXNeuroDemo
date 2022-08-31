%% testARMAXNeuro demonstrates a single parameter recovery run of the model. 
%  We generate simulated neural data according to certain parameters, then 
%  fit the model to determine if we can recover those parameters. Note that
%  there is a slight mismatch in the model used for simulation and the
%  model used for fitting because we don't know the value of 
%  some of the mean terms in the model prior to simulating. 
%  
%  If you change the simulation parameters, you may also have to change the
%  mean parameter. 
%  
%  Change which components are included in the model using the '_on' flags,
%  e.g., changing cue_mem_on to 0 would remove a cue memory component from
%  the model. 
%  
%  Change model parameters by changing tau/amp/coeff variables below. If
%  changing the model parameters results in a huge change in the 'mu'
%  value printed to the command line, then they should be changed back
%  because you may have changed them to unstable values, e.g., if you set
%  intrinsic_coeff(1) to 2, then there will be exponential growth in firing
%  rate across bins. 
%
%  One useful postprocessing step to improve model recovery could be to
%  only include components with high enough output.stat.delta_r2 value.
%  Essentially, only include components that explain enough variance. 

%% first, set simulation parameters
intrinsic_coeff = [.2, .1, .05, .05, 0]/5;
seasonal_coeff = [.2, .05, .05, .01, 0]/5;

cue_amp = -.5/2;
cue_tau = 10;

sample_amp = .8/2;
sample_tau = 20;

sample_exo = 2;
cue_exo = 2;

% which components to include in the moel
cue_mem_on = 1;
sample_mem_on = 1;
intrinsic_on = 1;
seasonal_on = 1;
exogenous_on = 1;

%% second, generaste random trial data
num_trials = 200;
sample = 2*(rand(num_trials,1)<.5)-1;
cue = 2*(rand(num_trials,1)<.5)-1;

%% third, generate spikes according to model
trial_length = 200; % trial length is 200 50 ms bins ms
trial_length_s = trial_length*50/1000; % trial length in seconds 
mean_fr = 10; % 10 hz
spike_prob = double(rand(num_trials,trial_length)<(mean_fr*50/1000));

% zero out all components
intrinsic_comp = 0;
seasonal_comp = 0;
cue_mem_comp = 0;
sample_mem_comp = 0;
cue_exo_comp = 0;
sample_exo_comp = 0;

% set predicted fr to initial fr
pred_fr = double(spike_prob);

cue_onset_t = 2;
sample_onset_t = 4;

% estimate mean parameter (may be dependent on choice of simulation parameters, so
% if you change those, you may have to change mu)
mu = 40;
for i = 1:num_trials
    for j = 1:trial_length
        if j>5
            if intrinsic_on
            intrinsic_comp = sum(spike_prob(i,j-5:j-1).*flip(intrinsic_coeff));
            end
        end
        if i>5
            if seasonal_on
            seasonal_comp = sum(spike_prob(i-5:i-1,j)'.*flip(seasonal_coeff));
            end
            if cue_mem_on
            cue_mem_comp = mu*sum(cue_amp*exp(-(trial_length_s*[1:5]+(.05*(j-1)))'/cue_tau).*flip(cue(i-5:i-1)));
            end
            if sample_mem_on
            sample_mem_comp = mu*sample_amp*sum(exp(-(trial_length_s*[1:5]+(.05*(j)))'/sample_tau).*flip(sample(i-5:i-1)));
            end
        end
        if exogenous_on
        if j*50/1000 > cue_onset_t
            cue_exo_comp = cue_exo*cue(i);
        else
            cue_exo_comp = 0;
        end
        if j*50/1000 > sample_onset_t
            sample_exo_comp = sample_exo*sample(i);
        else
            sample_exo_comp = 0;
        end
        end
        
        pred_fr(i,j) = mu+(nansum([intrinsic_comp,seasonal_comp, cue_mem_comp, sample_mem_comp, cue_exo_comp, sample_exo_comp]));
       
        p = pred_fr(i,j)*50/1000;
        o = 0;
        while p>rand()
            p = p-1;
            o = o+1;
        end
        %if j>5 | i > 1
        spike_prob(i,j) = o*1000/50;
        %mu = mean(spike_prob, 'all');
        %end
    end
end

disp("mu for next sim: " + mean(spike_prob, 'all'));

% translate spike rate into spike times
frs = spike_prob;
spikeTime = [];
signalTimeStart = [];
signalTimeEnd = [];
cueTimes = [];
sampleTimes = [];
for i = 1:num_trials
    for j =1:trial_length
        ct = trial_length_s*(i-1 + j/trial_length) + .01;
        p = frs(i,j)*50/1000;
        o = 0;
        while p>rand()
            p = p-1;
            o = o+1;
        end
        currTime = repmat(ct, [o,1]);
        spikeTime = [spikeTime; currTime];
    end
    signalTimeStart = [signalTimeStart; trial_length_s*(i-1)];
    signalTimeEnd = [signalTimeEnd; trial_length_s*i-.001];
    cueTimes = [cueTimes; trial_length_s*(i-1) + cue_onset_t + eps]; % cue onset 2 seconds after start of trial
    sampleTimes = [sampleTimes; trial_length_s*(i-1) + sample_onset_t + eps]; % sample onset 4 seconds after start of trial
end

%% fourth, construct neuron data input to the model
neuronData = NeuronData(spikeTime*1000, 50);
neuronData.addTrialBeginEndSignal('begin_trial', 'end_trial');
neuronData.addAlignTrialSignal('begin_trial', 100, 20000);
neuronData.addSignalTime('begin_trial', signalTimeStart*1000);
neuronData.addSignalTime('end_trial', signalTimeEnd*1000);
neuronData.addSignalTime('cue', cueTimes*1000);
neuronData.addSignalTime('sample', sampleTimes*1000);
neuronData.addSignalValue('cue', cue);
neuronData.addSignalValue('sample', sample);

%% fifth, construct model
ARMAXNeuroModel = ARMAXNeuro(neuronData, 0);

ARMAXNeuroModel.timeMaskSignal('cue_to_end', 'cue', 'end_trial');
ARMAXNeuroModel.timeMaskSignal('sample_to_end', 'sample', 'end_trial');

% Adding a short memory component with the depth of 5
short_interval = neuronData.binSize;
ARMAXNeuroModel.addMemShort(5, short_interval);

% Adding a seasonal memory component with the depth of 5
seasonal_interval = 10; % time between trials to use for computing seasonal timescales, 10 secs
ARMAXNeuroModel.addMemSeason(5, seasonal_interval);

% parameter limits for exogenous memory component
tau_upper = 60;
tau_lower = .5;
amp_lower = -4;
amp_upper = 4;
logn_params = [4, 1.67]; % parameters are drawn from a lognormal distribution

ARMAXNeuroModel.addMemExo('sample_mem', 'sample', 0, 'sample', 5, [amp_lower, amp_upper], [tau_lower, tau_upper], logn_params);    
ARMAXNeuroModel.addMemExo('cue_mem', 'cue', 0, 'cue', 5, [amp_lower, amp_upper], [tau_lower, tau_upper], logn_params);   

ARMAXNeuroModel.addExoSignal('sample', 'sample', 0, 'sample_to_end');
ARMAXNeuroModel.addExoSignal('cue', 'cue', 0, 'cue_to_end');

%% sixth, fit the model and do model selection
% Initializing the data model
ARMAXNeuroModel.initData();

crossIter = 5; % Cross-validation instances 
num_models = 32;
% cell array to save fitting results
[fitting_results, model_selection, output] = ARMAXNeuroModel.fitAllModelSelection(num_models, crossIter);

%% seventh, check if the modeling correctly recovered parameters
intrinsic_coeff = [.2, .1, .05, .05, 0]/5;
seasonal_coeff = [.2, .05, .05, .01, 0]/5;

cue_amp = -.5/2;
cue_tau = 10;

sample_amp = .8/2;
sample_tau = 20;

sample_exo = 2;
cue_exo = 2;

real.amp_cue_mem = cue_amp;
real.amp_sample_mem = sample_amp;
real.tau_cue_mem = cue_tau;
real.tau_sample_mem = sample_tau;
real.exo_sample = sample_exo;
real.exo_cue = cue_exo;
real.tau_intrinsic = ARMAXNeuroModel.getARTau(intrinsic_coeff, 50);
real.tau_seasonal = ARMAXNeuroModel.getARTau(seasonal_coeff, 10);

fields = fieldnames(real);
disp("---------------------------------")
disp("Model Recovery Results: ");
disp("---------------------------------")
disp("1. Does the model correctly recover the components used for simulating?" + ...
     " If so, the next two output lines should match.");
comps = ["exogenous", "intrinsic", "seasonal", "sample_mem", "cue_mem"];
compon = logical([exogenous_on, intrinsic_on, seasonal_on, sample_mem_on, cue_mem_on]);

% real model components
comp = comps(compon);
str = "model components- real: ";
for i = 1:length(comp)
    str = str + comp(i) + ", ";
end

disp(str);
% predicted components
comp = output.comp.labels(logical(output.comp.included));
str = "model components- pred: " ;
for i = 2:length(comp)
    str = str + comp(i) + ", ";
end
disp(str);

% parameters
disp("---------------------------------")
disp("2. Does the model correctly recover the parameters used for simulating?");
disp(" If estimate = NaN, then the component was not included in the model,");
disp(" so the predicted and real values need not match.");
disp(" Note that there is a slight mismatch in the model used for simulation ");
disp(" and the model used for fitting because we don't know the value of some ");
disp(" of the mean terms in the model prior to simulating. This may result in a ");
disp(" slight mismatch between the real/estimated values")

for i = 1:length(fields)
    disp(fields{i} + "- real: " + real.(fields{i}) + ", estimated " + num2str(output.var.(fields{i})));
end