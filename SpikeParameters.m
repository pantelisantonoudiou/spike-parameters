classdef SpikeParameters < matlab.mixin.Copyable
    %
    % Start object instance
    % x = SpikeParameters(store_mat{10,1}(1,:), 10000)
    %
    % Get spike parameters
    % properties = x.spike_parameters(false)
    
    properties
        interp_factor = 10                  % interpolation factor
        data_trim_ms = 25                   % milliseconds
        ap_dt_threshold = 10                % action potential gradient threshold
        
        % phase plane
        trim_spikes = [1, 1]                % spikes to remove from phase_plane [start, stop]
        spike_window = [3, 3]               % time to extract from spike peak [before, after] (msec)
        spike_window_samp                   % time to extract from spike peak [before, after] (samples)
        
        % calculated values
        data                                % data
        t                                   % time vector
        Fs                                  % sampling rate (per second)
        data_trim                           % data points
        min_peak_dist                       % minimum peak distance points
    end
    
    methods(Static)
        % plot phase plane
        function p = get_phase_plot(spikes, mean_col, fill_col)
            % p = get_phase_plot(spikes, mean_col, fill_col)
            
            % get average action potential
            avg_spike = mean(spikes,1);
            
            % get gradient
            gradient_spikes = gradient(spikes);
            
            % find average and standard deviation
            avg_gradient_spikes = mean(gradient_spikes,1);
            sem_gradient_spikes = std(gradient_spikes,1)/sqrt(size(spikes,1));
            
            % get error
            error_pos = avg_gradient_spikes + sem_gradient_spikes;
            error_neg = avg_gradient_spikes - sem_gradient_spikes;
            
            % get x and y for plotting
            x = avg_spike;
            y = avg_gradient_spikes;
            
            % plot
            hold on
%             xfill = horzcat(x, fliplr(x));
%             yfill = horzcat(error_pos, fliplr(error_neg));
%             fill(xfill, yfill, fill_col,'LineStyle','none','DisplayName','SEM');
            p = plot(x, y,'color', mean_col, 'visible', 'off', 'linewidth', 1.5);
            plot(reshape(spikes',[],1), reshape(gradient_spikes',[],1),'color', mean_col);
        end
        
        % plot average waveform
        function p = aver_spike_waveform(spikes, Fs,  mean_col, fill_col)
            % p = aver_spike_waveform(spikes, Fs,  mean_col, fill_col)
            % spikes = matrix where rows = different spikes and
            % cols = time
            
            % get average action potential
            aver_spike = mean(spikes,1);
            
            % get standard deviation of action potentials
            sem_spike = std(spikes,1)/sqrt(size(spikes,1));
            
            % get negative and positive error
            y_error_neg = aver_spike - sem_spike;
            y_error_pos = aver_spike + sem_spike;
            
            % get x-axis time
            x = (0:size(spikes, 2)-1)/Fs*1000;
            y = aver_spike;
            
            % plot mean and shaded sem
            hold on;
%             xfill = horzcat(x, fliplr(x));
%             yfill = horzcat(y_error_pos, fliplr(y_error_neg));
%             fill(xfill, yfill, fill_col,'LineStyle','none','DisplayName','SEM');
            p = plot(x, y,'color', mean_col,'Linewidth',1.5, 'visible', 'off');
            plot(x,spikes,'color', mean_col)
        end
    end
    
    methods
        
        % consutructor (interpolate and get properties)
        function obj = SpikeParameters(data, Fs)
            
            % get interpolation sampling rate
            obj.Fs= Fs*obj.interp_factor;
            
            % create interpolation x vectors
            t = (0:1:length(data)-1)/obj.Fs';
            obj.t = (0:1/obj.interp_factor:length(data)-1)/obj.Fs';
            
            % interpolate data
            obj.data = interp1(t, data, obj.t, 'spline');
            
            % calculate number of points to trim
            obj.data_trim = obj.data_trim_ms/1000 * obj.Fs;
            
            % set min peak distance to 1.5 ms
            obj.min_peak_dist = round(0.0015 * obj.Fs);
            
            % get time to extract spikes in samples
            obj.spike_window_samp = round(obj.spike_window * obj.Fs / 1000);
        end
        
        % Get spike locations
        function [pks, locs] = spike_detect(obj, threshold)
            % [pks, locs] = obj.spike_detect()
            
            if nargin == 1
                % spike threshold (gradient crossing speed)
                % threshold = 15; %20 (30 too stingent for rheobase)
                threshold = 5*std(gradient(obj.data(obj.data_trim:end-obj.data_trim))*obj.Fs);
                if threshold < 5
                    threshold = 5;
                elseif threshold > 20
                    threshold = 20;
                end
            end
            
            % detect rising peaks
            idx = peakseek(gradient(obj.data)*obj.Fs,obj.min_peak_dist,threshold);
            
            if isempty(idx)==1 % if no peaks detected exit
                pks = [];
                locs = [];
                return
            end
            
            % preallocate vector
            locs = zeros(1,length(idx));
            
            for i = 1:length(idx)
                [~,idx2] = max(obj.data(idx(i):idx(i) + obj.min_peak_dist));
                locs(i) = idx(i)+ idx2 - 1;
            end
            
            % get peaks
            pks = obj.data(locs);
            
            % remove peaks smaller than the median value
            locs = locs(pks > median(obj.data));
            pks = pks(pks > median(obj.data))';
            
            %              plot(obj.data);hold on
            %              plot(locs, pks, 'rx')
        end
        
        % Get spike properties
        function properties = spike_parameters(obj, plot_var)
            % export spike parameters
            
            % Get peak locations and amplitude
            [pks, locs] = obj.spike_detect();
            
            % Preallocate storage vectors
            APthresh_idx = zeros(size(locs));
            AHP_idx = zeros(size(locs)); post = obj.Fs*0.01;
            for i = 1:length(locs)
                
                % Find AP Threshold
                start_point = locs(i) - obj.min_peak_dist; % get template start point
                spike_template = obj.data(start_point:locs(i));
                APthresh_idx(i) = start_point + find(gradient(spike_template)*obj.Fs> obj.ap_dt_threshold, 1, 'first');
                
                % Find AHP
                [~, idx] = min(obj.data(locs(i):locs(i)+post));
                AHP_idx(i) = locs(i) + idx;
            end
            
            % Find AP amp, AHP amp
            AP_amp = pks' - obj.data(APthresh_idx);
            AHP_amp = obj.data(APthresh_idx) - obj.data(AHP_idx);
            
            % Find Half width
            AP_wdth_idx = zeros(2, length(locs));
            for  i = 1:length(locs)
                % Get half amplitude
                half_amp = obj.data(locs(i)) - (abs(AP_amp(i))/2);
                [~,idx1] = min(abs(obj.data(APthresh_idx(i):locs(i))- half_amp)); % find half amp index left of peak
                [~,idx2] = min(abs(obj.data(locs(i):locs(i) + obj.min_peak_dist)- half_amp)); % % find half amp index right of peak
                
                % Get Half width index
                AP_wdth_idx(1,i) = APthresh_idx(i) + idx1 + 1;
                AP_wdth_idx(2,i) = locs(i) + idx2 - 1 ;
            end
            
            %%% ------------------------- Plot ------------------------- %%%
            if plot_var == true
                plot(obj.t, obj.data,'k'); hold on;
                plot(obj.t(locs),obj.data(locs),'rx','LineWidth',2)
                plot(obj.t(AP_wdth_idx),obj.data(AP_wdth_idx),'bx','LineWidth',2)
                plot(obj.t(APthresh_idx),obj.data(APthresh_idx),'gx','LineWidth',2)
                plot(obj.t(AHP_idx), obj.data(AHP_idx),'mx','LineWidth',2)
                set(gca,'Box','off','FontName', 'Arial','FontSize', 12);
                set(gca,'TickDir','out'); set(gcf,'color','w');
            end
            %%% -------------------------------------------------------- %%%
            
            %%% ---------------- Get values for export ----------------- %%%
            properties.AP_amp =  AP_amp * 1000;                                                 % AP amplitude (mV)
            properties.AHP_amp  =  AHP_amp * 1000;                                              % AHP amplitude (mV)
            properties.AP_threshold = obj.data(APthresh_idx) * 1000;                            % AP threshold (mV)
            properties.isi = diff(locs)/obj.Fs*1000;                                            % inter spike interval (ms)
            properties.peak_to_trough = (AHP_idx - locs) / obj.Fs * 1000;                       % Peak to trough (ms)
            properties.AP_half_width = (AP_wdth_idx(2,:)- AP_wdth_idx(1,:)) / obj.Fs * 1000;    % Half width (ms)
            %%% -------------------------------------------------------- %%%
        end
        
        % Extract spikes from data
        function spikes = extract_spikes(obj)
            
            % Get peak locations and amplitude
            [~, locs] = obj.spike_detect();
            
            % discard spikes at onset and offset of current stimulation
            locs = locs(1 + obj.trim_spikes(1) : end - obj.trim_spikes(2));
            
            % create matrix to store spikes
            spikes = zeros(length(locs), obj.spike_window_samp(1) + obj.spike_window_samp(2));
            
            for i = 1:length(locs)
                spikes(i,:) = obj.data(locs(i)-obj.spike_window_samp(1) : locs(i)+obj.spike_window_samp(2)-1);
            end
            
            
        end
  
    end
    
end






