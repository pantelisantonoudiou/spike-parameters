% User input for phase plane and spike average waveforms
% -------------------------
data_path = 'Z:\Phil\Ephys\phaseplane_test';
spike_array_path = 'Z:\Phil\Ephys\phaseplane_test\spike_array.mat';
filters = {'PV_Sa'}; %{'PV_K', 'PV_Sa'}
Fs = 10000;
plot_var = 'mean'; % 'mean' or 'ind'
% -------------------------

% load spike_array
s = load(spike_array_path);

% init object
spike_obj = AllCells(s.spike_array, Fs, plot_var);

% 1- Get phase plane plot
spike_obj.plot_phase_plane_conditions(data_path, filters)

% 2- Get waveform
spike_obj.plot_spike_waveform_conditions(data_path, filters)
    