% User input for phase plane and spike average waveforms
% -------------------------
path = 'Z:\Phil\Ephys\phaseplane_test';
filters = {'PV_K', 'PV_Sa'};
Fs = 10000;
plot_var = 'mean'; % 'mean' or 'ind'
% -------------------------

% init object
spike_obj = AllCells();
spike_obj.Fs = Fs;
spike_obj.plot_var = plot_var;

% 1- Get phase plane plot
spike_obj.plot_phase_plane_conditions(path, filters)

% 2- Get waveform
spike_obj.plot_spike_waveform_conditions(path, filters)
    