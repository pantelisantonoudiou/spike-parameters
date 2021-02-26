classdef AllCells < matlab.mixin.Copyable
    %
    % prism_array = AllCells.all_cells('Z:\Pantelis\Phil_data')
    %
    
    properties
    end
    
    methods(Static)
        
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
end






