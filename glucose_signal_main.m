% ===================================================================== %  
% ============= TUE009 Glucose sensor data preprocessing ============== %
% ===================================================================== % 

% Author: María Berjano, 2023
%
% This script processes the data for multiple participants and generates 
% various outputs, including adjusted signals, summary matrices, binned 
% data, and a quality control (QC) table. 
%
% READ IN FILES
% Subjects glucose data. The folder structure needs to be the following:
% - current folder/P0XX/GlucoseP0XX_glucose_DD-M-YYYY.csv
% where "current folder" is the directory where glucose_signal_main.m is
% stored, and XX is the subject ID
% e.g."...\Glucose_Sensor\P038\GlucoseP038_glucose_16-12-2022.csv"
%
% REQUIRES:	
% The required program files for this script are: 
% - MATLAB 2020a or later
% - Bioinformatics Toolbox
% The required functions for this script are:
% - inferno.m
% - msbackadj_mod.m
% OUTPUT
% TG_all_subjects.csv

close all; clear; clc
scriptsDir = pwd;
parentDir = fullfile(pwd, '..');
cd(parentDir)

% Get list of participants' IDs
listing = dir();
listing = listing(~ismember({listing.name}, {'.', '..'})); % Exclude the current directory (.) and parent directory (..)
IDs = char({listing([listing.isdir] & contains({listing.name}, 'P0') & cellfun('length',{listing.name})==4).name}'); % Filter the listing (exclusion of subject 87 is temporary)
path = pwd;
init_path = extractBefore(path,'Projects');

% Define parameters to bin the data in uniform timesteps
fs = minutes(5); % sampling frequency
window_length = hours(1); 
window_stride = minutes(5);

% Define time range to look for the baseline
start_time = 0; 
end_time = 8;

% Loop through participants list
for i = 1:length(IDs)
    
    % ============ Create timetable with the recorded values ============
    id = IDs(i,:);

    % Access table with glucose values
    T = dir([id, '/*.csv']);
    idx = ~cellfun('isempty', regexp(lower({T.name}), 'glucose')); % Perform case-insensitive matching
    T = T(idx);
    if length(T) > 1
        
        fileDates = [T(:).datenum]; % Get the dates of the matching files
        [~, maxIdx] = max(fileDates); % Find the index of the file with the latest date
        T = T(maxIdx,:); % Retrieve the file with the latest date

    end
    T = char(T.name);
    data = readtable([path,filesep,id,filesep,T]);

    % Take only the data with recording type 0
    data = data(data.Aufzeichnungstyp == 0,:);

    % Convert glucose units if necessary
    Glucose = data.Glukosewert_VerlaufMg_dL;
    if any(Glucose < 10)
        Glucose = Glucose * 18.02; % mmol/l to mg/dl
    end

    Time = datetime(data.Ger_tezeitstempel,'InputFormat','dd-MM-yyyy HH:mm','TimeZone','local');
    TG = sortrows(timetable(Time,Glucose));
    TG = TG(year(TG.Time)>=2022,:);

    % Remove duplicate recordings (two rows with the same date and time) if
    % necessary
    if length(unique(TG.Time)) < height(TG) % Check if there are duplicate times
        
        uniqueTimes = unique(TG.Time); % Find unique times
    
        % Preallocate a table to store averaged values
        TG_unique = table(uniqueTimes, zeros(size(uniqueTimes)), 'VariableNames', {'Time', 'Glucose'});
    
        for j = 1:length(uniqueTimes) % Loop through each unique time and calculate the average Glucose
            TG_unique.Glucose(j) = mean(TG.Glucose(TG.Time == uniqueTimes(j)));
        end

        TG = TG_unique;
    end % If there are no duplicates, just use the original table

    % ========================== Get estimated baseline ===========================
    [baseline_value(i), baseline_time{i}] = Get_baseline_offset(TG, start_time, end_time, fs, window_length, window_stride);

    % ================ Get adjusted and smoothed signal =================
    addpath (parentDir)
    addpath (scriptsDir)
    [full_signal,TG_adj] = Get_adjusted_signal(TG,id,baseline_value(i));   

    % Add estimated intravascular glucose value
    TG_adj.EstIVGlucose = interp1(TG_adj.Time, TG_adj.Glucose, TG_adj.Time + minutes(5), 'linear', 'extrap');

    % =============== Get binned data and summary matrix ================
    [matrix, mat_fig, total_Days, TG_binned] = Get_matrix(TG_adj, window_length, window_stride, id, baseline_value(i), baseline_time{1,i});

    % =========================== Save output ===========================
    % Plots
    
    saveas(full_signal,[path,filesep,id,filesep,strcat(id,'_signal_adjusted.png')])
    saveas(full_signal,[init_path,'Projects\TUE009_DFG_Glucose_RPE\11_Interim_analyses\Plot_Glucose_Sensor\Summary_plots\',strcat(id,'_signal_adjusted.png')])
    saveas(mat_fig,[path,filesep,id,filesep,strcat(id,'_summary_matrix_adjusted.png')])
    saveas(mat_fig,[init_path,'Projects\TUE009_DFG_Glucose_RPE\11_Interim_analyses\Plot_Glucose_Sensor\Summary_plots\',strcat(id,'_summary_matrix_adjusted.png')])
    % Matrices
    save([path,filesep,id,filesep,strcat(id,'_signal_adjusted.mat')],"TG_adj")
    save([path,filesep,id,filesep,strcat(id,'_summary_adjusted.mat')],"matrix")
    % Timetable
    save([path,filesep,id,filesep,strcat(id,'_binned_data.mat')],"TG_binned")
    
    
    % ==================== Compute parameters for QC ====================
    MissingDataPercent(i,1) = 100*nnz(isnan(matrix.glucose))/(144*height(matrix.glucose));
    MissingDataDays(i,1) = (10*nnz(isnan(matrix.glucose)))/(60*24);
    DaysRecorded(i,1) = length(total_Days);
    MinVal(i,1) = min(matrix.glucose,[], 'all');
    MaxVal(i,1) = max(matrix.glucose,[], 'all');

    close all

end

% ======== Create table with glucose recordings from all subjects ========
BaselineValue = baseline_value';
BaselineTime = baseline_time';
Summary_Table = table(IDs, DaysRecorded, MissingDataPercent, MissingDataDays, MinVal, MaxVal, BaselineValue, BaselineTime);
[B,Outlier] = rmoutliers(Summary_Table.DaysRecorded); % outlier: more than three scaled median absolute deviations (MAD) 
Summary_Table = addvars(Summary_Table,Outlier,'Before',"MissingDataPercent");
writetable(Summary_Table,[init_path,'Projects/TUE009_DFG_Glucose_RPE/11_Interim_analyses/Plot_Glucose_Sensor/Summary_plots/qc_table.csv'])

% Save table with all glucose and timestamps of all participants
for i = 1:length(IDs)

    % Get glucose data from participant and create timetable with the
    % recorded values
    id = IDs(i,:);
    load(strcat(path,filesep,id,filesep,id,"_binned_data.mat"));
    subj_table{i} = cell2table(cell(height(TG_binned),4)); %Generate table with required rows and columns
    subj_table{i}.Properties.VariableNames = {'id','time','gl','estIVglu'}; %Give variable names
%     subj_table{i}.id(:) = {strcat('Subject',{' '},string(str2double(extractAfter(id,'P0'))))};
    subj_table{i}.id = ones(height(TG_binned),1)*str2double(extractAfter(id,'P0'));
    subj_table{i}.time = TG_binned.Time;
    subj_table{i}.gl = TG_binned.Glucose;
    subj_table{i}.estIVglu = TG_binned.EstIVGlucose;

end

TG_all_subjects = vertcat(subj_table{:});
TG_all_subjects = rmmissing(TG_all_subjects);
writetable(TG_all_subjects,[strcat(path,filesep,"TG_all_subjects.csv")])

%% Functions

function [glucose_baseline, time_baseline] = Get_baseline_offset(TG, start_time, end_time, fs, window_length, window_stride)
        
    % Identifies the baseline glucose level of a glucose sensor recording
    %       
    % Parameters:
    %     TG: timetable
    %          time and glucose values
    %     start_time: double
    %          time of the day to start searching for baseline value 
    %     end_time: double
    %          time of the day to stop searching for baseline value   
    %     fs: duration
    %          sampling frequency
    %     window_length: duration
    %          length of the window used to search the baseline
    %     window_stride: duration
    %          stride of the window used to search the baseline
    % Returns: 
    %     glucose_baseline: double
    %          estimated glucose baseline value
    %     time_baseline: duration
    %          time of the day of the estimated glucose baseline level

    % Filter values at the specified day period
    TG_filtered = TG(hour(TG.Time) <= end_time & hour(TG.Time) >= start_time,:);
    
    reach_end = 0;
    n = 0;
    t_current = min(timeofday(TG_filtered.Time));
    n_samples_W = length(unique(TG_filtered.Time.Day))*(window_length/fs); % ideal number of samples per window
    
    % Slide window over the selected search time period
    while reach_end == 0
            
            % Get values of the current window
            TG_window = TG_filtered(timeofday(TG_filtered.Time) < (t_current + window_length) & timeofday(TG_filtered.Time) >= t_current,:);
            
            % compute parameters only if there are more than 80% the
            % samples in the window
            if (n_samples_W - height(TG_window))/n_samples_W < 0.2

                n = n + 1;
                B = rmoutliers(TG_window.Glucose); % remove outliers
                window_std(n) = std(B); % compute standard deviation
                window_mean(n) = mean(B); % compute average
                window_time(n) = t_current;

            end
            
            % slide window to the next position
            t_current = t_current + window_stride;
            
            % check if the window reached the end of the search time period
            if t_current + window_length >= max(timeofday(TG_filtered.Time))
                
                reach_end = 1;
            
            end
        
    end
    
    % Define the baseline as the average of the glucose values in the
    % window with lowest standard deviation
    glucose_baseline = min(window_mean(window_std == min(window_std)));
    time_baseline = min(window_time(window_mean == glucose_baseline));

end

function [full_signal,TG_adj] = Get_adjusted_signal(TG,id,baseline_val)
    
    % Corrects baseline of original glucose signal and adds an estimated
    % baseline value as offset
    %       
    % Parameters:
    %     TG: timetable
    %          time and glucose values
    %     id: char or string
    %          subject's identification number
    %     baseline_val: double
    %          estimated glucose baseline value 
    % Returns: 
    %     full_signal: figure
    %          figure of complete adjusted signal and baseline
    %     TG_adj: timetable
    %          adjusted timetable

    % =============== Correct baseline of original signal =============== 
    % Define different sections in the recording. Consider different
    % sessions if there is a gap of more than 20 minutes in the recording
    sections = [0,find(diff(TG.Time)>minutes(20))',numel(TG.Time)];
    adj_signal = nan([length(TG.Glucose),1]);
    baseline_reg = nan([length(TG.Glucose),1]);

    % Correct the signal separately for each section using a slightly
    % modified version of the msbackadj function
    for i = 2:numel(sections)

        glucose_sections = TG.Glucose(sections(i-1)+1:sections(i));
        [b, yOut] = msbackadj_mod((1:length(glucose_sections))',glucose_sections,'WindowSize',300,'StepSize',100,'ShowPlot',false);%'PreserveHeights',true);
        adj_signal(sections(i-1)+1:sections(i)) = yOut;
        baseline_reg(sections(i-1)+1:sections(i)) = b;

    end
    
    % ============= Add offset to corrected glucose signal ==============
    TG_adj = sortrows(timetable(TG.Time,(adj_signal + baseline_val)));
    TG_adj = renamevars(TG_adj,'Var1','Glucose');
    
    % ====================== Smooth glucose signal ======================
    TG_adj.Glucose = smoothdata(TG_adj.Glucose,'gaussian',5);

    % Plot
    full_signal = figure;
    hold on
    for i = 2:numel(sections)
        plot(TG_adj.Time(sections(i-1)+1:sections(i)), TG_adj.Glucose(sections(i-1)+1:sections(i)), 'k')
    end
    title(strcat("Glucose level summary ",id))
    subtitle(strcat("Baseline glucose level = ", num2str(round(baseline_val)), " mg/dL"),'FontSize',8)
    ylabel("Glucose level (mg/dL)")
    yline(baseline_val,'m')
    hold off
    
end

function [matrix, mat_fig, total_Days, TG_binned] = Get_matrix(TG, window_length, window_stride, id, baseline_val, baseline_t)
    
    % Generates a matrix from a glucose recording in which the rows
    % correspond to the days and the columns to the time bins in a day
    %       
    % Parameters:
    %     TG: timetable
    %          time and glucose values
    %     window_length: duration
    %          length of the window used to search the baseline
    %     window_stride: duration
    %          stride of the window used to search the baseline
    %     id: char or string
    %          subject's identification number
    %     baseline_val: double
    %          estimated glucose baseline value 
    %     baseline_t: duration
    %          time of the day of the estimated glucose baseline level
    % Returns: 
    %     matrix: struct
    %          includes the glucose matrix (double) and general information
    %          about the subjects's recording (table)
    %     mat_fig: figure
    %          figure of resulting matrix
    %     total_Days: datetime
    %          number of recorded days
    %     TG_binned: timetable
    %          binned time and glucose values

    % Define bins and initialize matrix - The matrix has as many rows as
    % days recorded and as many columns as the number of bins that fit in a
    % day taking into account the window stride
    firstDate = dateshift(groupsummary(TG.Time,findgroups(year(TG.Time)),@min),'start', 'day');
    lastDate = dateshift(groupsummary(TG.Time,findgroups(year(TG.Time)),@max),'start', 'day');
    total_Days = datetime(firstDate:lastDate,'Format','dd-MM-yyyy');
    bins = (timeofday(datetime(0,0,0,'Format','HH:mm:ss')) + minutes(0:minutes(window_stride):24*60-1));
    matrix.glucose = NaN([length(total_Days),length(bins)]);
    matrix.info = array2table([string(datetime(firstDate,'Format','dd-MM-yyyy')),string(datetime(lastDate,'Format','dd-MM-yyyy')),string(minutes(window_stride)),baseline_val, string(baseline_t),string(hours(window_length)), string(minutes(window_stride))],"VariableNames",["First Day","Last Day","Bin Size (min)","Baseline_val","Baseline_time","Window_length","Window_stride"]);
    matrixIV.glucose = NaN([length(total_Days),length(bins)]);
    matrixIV.info = array2table([string(datetime(firstDate,'Format','dd-MM-yyyy')),string(datetime(lastDate,'Format','dd-MM-yyyy')),string(minutes(window_stride)),baseline_val, string(baseline_t),string(hours(window_length)), string(minutes(window_stride))],"VariableNames",["First Day","Last Day","Bin Size (min)","Baseline_val","Baseline_time","Window_length","Window_stride"]);
    
    % Fill matrix with the glucose level correspponding to each cell (bin)
    for d = 1:length(total_Days)
        data_d = TG(day(TG.Time,'dayofyear') == day(total_Days(d),'dayofyear'),:);
        for b = 1:length(bins)
            data_bin = data_d(timeofday(data_d.Time) < bins(b) + window_stride & timeofday(data_d.Time) >= bins(b),:);
            if ~isempty(data_bin)
                matrix.glucose(d,b) = mean(data_bin.Glucose);
                matrixIV.glucose(d,b) = mean(data_bin.EstIVGlucose);
            end
        end
    end

    % Plot baseline on top of the summary matrix
    mat_fig = figure;
    a = imagesc(matrix.glucose);
    % 1:12:length(bins) before instead of  1:(120/minutes(window_stride)):length(bins)
    set(gca,'xtick', 1:(120/minutes(window_stride)):length(bins),'xticklabel',string(0:2:23),'FontSize',10)
    set(gca,'ytick',1:5:length(total_Days),'yticklabel',string(1:5:length(total_Days)),'FontSize',10)
    xlabel("Time (bin: 5 min)")
    ylabel("Day")
    title(strcat("Glucose level summary ",id),'FontSize',10)
    subtitle(strcat("Baseline glucose level = ", num2str(round(baseline_val)), " mg/dL"),'FontSize',8)
    set(a,'AlphaData',~isnan(matrix.glucose))
    colormap(inferno)
    colorbar
    hold on
    x = baseline_t-bins;
    [~,x1] = min(x(x>=0));
    rectangle('Position',[x1-0.5,0.5,(window_length/window_stride),length(total_Days)],'LineStyle','-','EdgeColor','m')
    hold off

    % Save timetable with starting time and glucose value of each bin
    total_bins = [];
    for d = 1:length(total_Days)
        total_bins = [total_bins datetime(total_Days(d)+bins,'Format','dd-MM-yyyy HH:mm:ss')];
    end
    TG_binned = sortrows(timetable(total_bins',reshape(matrix.glucose.',1,[])',reshape(matrixIV.glucose.',1,[])'));
    TG_binned = renamevars(TG_binned,["Var1","Var2"],["Glucose","EstIVGlucose"]);
    
end