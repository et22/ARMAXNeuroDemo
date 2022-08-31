# ARMAXNeuro
## Input
The input to ARMAXNeuro is a `NeuronData` object that contains signal times and values. The construction of the NeuronData object is demonstrated in testARMAXNeuro. To construct the object, first call the constructor with input `spikeTime`, containing spike times in milliseconds and `binSize` containing length of bin in milliseconds for model. 
```
neuronData = NeuronData(spikeTime*1000, 50);
```
Next, call `addTrialBeginEndSignal` to label the signal times associated with the beginning and end of the trial.
```
neuronData.addTrialBeginEndSignal('begin_trial', 'end_trial');
```
Then, call `addAlignTrialSignal` to label the signal used for aligning time bins, e.g., for seasonal timescales. Also provide the maximum time from the trial align signal to the end of the previous trial (here, 100) and the maximum time from the trial align signal to the end of the current trial (here, 20000). 
```
neuronData.addAlignTrialSignal('begin_trial', 100, 20000);
```
Finally, add signal times (in ms) and values to neuron data. An example with adding cue and sample times and values is shown below. 
```
neuronData.addSignalTime('cue', cueTimes*1000);
neuronData.addSignalTime('sample', sampleTimes*1000);
neuronData.addSignalValue('cue', cue);
neuronData.addSignalValue('sample', sample);
```
## Output
The `fitAllModelSelection` method returns three outputs - [fitting_results, model_selection, output]. The output structure should contain all of the information you need from fitting. The output structure has the following fields:  
- **var** - contains the model parameter estimates 
- **p** - contains corresponding p values
- **stat** - contains R2 of best model, and R2 associated with removing each component from the full model
- **comp** - lists which components were included in the best fitting model
## Testing
The output from an example run of testARMAXNeuro are shown below. They should be self explanatory. 
```
---------------------------------
Model Recovery Results: 
---------------------------------
1. Does the model correctly recover the components used for simulating? If so, the next two output lines should match.
model components- real: exogenous, intrinsic, seasonal, sample_mem, cue_mem, 
model components- pred: exogenous, intrinsic, seasonal, sample_mem, cue_mem, 
---------------------------------
2. Does the model correctly recover the parameters used for simulating?
 If estimate = NaN, then the component was not included in the model, so the predicted and real values need not match.
 Note that there is a slight mismatch in the model used for simulation 
 and the model used for fitting because we don't know the value of some 
 of the mean terms in the model prior to simulating. This may result in a 
 slight mismatch between the real/estimated values
amp_cue_mem- real: -0.25, estimated -0.14976
amp_sample_mem- real: 0.4, estimated 0.27002
tau_cue_mem- real: 10, estimated 11.6189
tau_sample_mem- real: 20, estimated 21.6388
exo_sample- real: 2, estimated 2.0203
exo_cue- real: 2, estimated 2.0272
tau_intrinsic- real: 50.0526, estimated 58.4806
tau_seasonal- real: 7.9634, estimated 9.4773
```