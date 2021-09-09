classdef AllCells < matlab.mixin.Copyable
    %
    % prism_array = AllCells.all_cells('Z:\Pantelis\Phil_data')
    %
    
    properties
        Fs = 10000
    end
    
    methods(Static) % Spike properties
        
        % get all spike properties in prism array
        function prism_array = all_cells(main_path)
            % set variables
            Fs = 10000;
            spike_train = 15;
            
            % get directory
            file_dir = dir(fullfile(main_path, '*mat'));
            
            % create storage array
            prism_array = cell(2,length(file_dir));
            prism_array(1,:) = {file_dir.name};
            
            % init progress bar
            w = waitbar(0, 'Please wait');
            for i = 1:length(file_dir)
                
                % load spike trains
                s = load(fullfile(main_path, file_dir(i).name),'store_mat');
                
                % get spike properties
                prism_array{2,i} = AllCells.repetitions(s.store_mat{spike_train,1}, Fs);
                
                waitbar(i/length(file_dir), w, 'Extracting Power...'); % update progress bar
            end
            
            close(w) % close progress bar
        end
        
        % get spike properties for one cell
        function neuron_array = repetitions(spike_matrix, Fs)
            % neuron_array = repetitions(spike_matrix, Fs)
            
            %%% ------------------ Get spike properties ------------------ %%%
            for i = 1:size(spike_matrix,1)  % iterate through repetitions
                
                % init SpikeParameters object and get properties
                x = SpikeParameters(spike_matrix(i,:), Fs);
                properties = x.spike_parameters(false);
                
                if i == 1 % create storage array
                    n_spike_properties = numel(fieldnames(properties));             % number of spike properties
                    temp_array = cell(size(spike_matrix,1), n_spike_properties);    % temporary array to store data
                    neuron_array = cell(2, n_spike_properties);                     % final array
                    neuron_array(1,:) = fieldnames(properties);                     % pass names of properties to array
                end
                
                % get data into cell array format
                temp_array(i,:) = struct2cell(properties);
            end
            
            %%% ---------- Restructure to format for PRISM --------------- %%%
            for ii = 1:n_spike_properties
                
                % create matrix for restorage
                repeat_matrix = NaN(500, size(spike_matrix, 1));
                for i = 1:size(spike_matrix, 1)
                    % pass data to new_matrix
                    repeat_matrix(1:length(temp_array{i,ii}), i) = temp_array{i,ii};
                end
                
                % remove rows where all nans
                neuron_array{2, ii} = repeat_matrix(any(~isnan(repeat_matrix),2),:);
            end
        end
        
    end
    
    methods(Static) % Phase plane
        
        function plot_conditions(main_path, conditions, Fs)
            
            figure()
            for i = 1:length(conditions)
                                
                % get spikes
                spikes = AllCells.aver_phase_plane(main_path, conditions{i}, Fs);
                
                % plot phase plane
                [col_mean,col_sem] = color_vec(i + 1); % get colors
                p(i) = SpikeParameters.get_phase_plot(spikes, col_mean, col_sem);
                
%                 p(i) = SpikeParameters.aver_spike_waveform(spikes, Fs, col_mean, col_sem);
                
            end
            
            prettify(gca)
%             xlabel('Time (ms)')
%             ylabel('Vm')
            
            xlabel('Vm (mV)')
            ylabel('dV/dt')
            legend(p, strrep(conditions, '_', ' '))

            
        end
        
        % get average spike for each cell
        function spikes = aver_phase_plane(main_path, filter, Fs)
            % spikes = aver_phase_plane(main_path, filter, Fs)
            %
            % INPUTS
            % ----------
            % main_path: str, path to folder with mat data
            % filter: str, used to filter files
            % Fs: int, sampling rate
            %
            % OUTPUTS
            % ----------
            % spikes: mat, rows = aver spike waveform form each cell
            % cols = time
            
            % init progress bar
            w = waitbar(0, 'Please wait');
            
            % get file directory
            file_dir = dir(fullfile(main_path, horzcat('*', filter, '*')));
            
            % create empty cell array
            all_spikes = cell(length(file_dir),1);
            for i = 1:length(file_dir)
                
                % load all current steps
                s = load(fullfile(main_path, file_dir(i).name ),'store_mat');
                
                % get index for largest current step
                io = cell2mat(s.store_mat(:,2));
                data = s.store_mat{get_index(s.store_mat(:,2), max(io)), 1};
                
                % get aver spike
                spikes = AllCells.choose_repetition(data, Fs, 0);
                all_spikes{i} = mean(spikes,1);
                
                waitbar(i/length(file_dir), w, 'Extracting Waveforms...'); % update progress bar
            end
            
            close(w) % close progress bar
            
            % combine all spikes
            spikes = vertcat(all_spikes{:});
            
        end
        
        % extract spike waveform
        function spikes = choose_repetition(data, Fs, plot_var, varargin)
        % spikes = choose_repetition(data, Fs, plot_var, varargin)
        % spikes = choose_repetition(data, 10000, 1, 2)
            
            % combine all sweeps
            if isempty(varargin) == 1
                all_spikes = cell(size(data, 1),1);
                for i = 1:size(data, 1) % iterate over sweeps
                    
                    % extract spikes
                    x = SpikeParameters(data(i,:), Fs);
                    all_spikes{i} = x.extract_spikes();
                end
                
                % combine all spikes
                spikes = vertcat(all_spikes{:});
            else
                % get sweep number from user input
                rep = varargin{1};
                
                % check if repetition number is within bounds
                if rep == 0 || rep > size(data,1)
                    error([' Repetition -' num2str(rep)  '- was not found.'])
                else
                    % extract spikes from sweep
                    x = SpikeParameters(data(rep,:), Fs);
                    spikes = x.extract_spikes();
                end
            end
            
            % plot spikes
            if plot_var == 1
                SpikeParameters.get_phase_plot(spikes)
            end
            
        end
        
    end
end






