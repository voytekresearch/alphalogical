%% Data Extraction for Erik
%

%% Imports

% Add eeglab to path and launch
addpath('/Users/thomasdonoghue/Documents/Software/eeglab13_4_4b');
eeglab;

%% Settings

% rtPB-3 Data Path
rtPB_data_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/Phase/Experiments/rtPB/2-Data/rtPB-3/processed/EEG/';

% PBA Data Path
pba_data_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/Phase/Experiments/PBA/2-Data/PBA-3/processed/EEG/';

% Path to save extracted data to
save_path = '/Users/thomasdonoghue/Desktop/ErikData/';

% Event code labels
rtPB_rest_ec = {'StartRest'};
pba_rest_ec = {'Rest_Start'};
rtPB_trial_ec = {'Start Block'};
pba_trial_ec = {'Exp_Block_Start'};

%% rtPB Data Extraction

% Get list of all available subject files
rtPB_files = get_files(rtPB_data_path, 'set');

% Set up for save name
num = 1000;

% Loop through all subjects, extracting and saving out data
for subj = rtPB_files
    
    s_name = [num2str(num), '.mat'];
    extract_save_data(subj, rtPB_data_path, rtPB_rest_ec, rtPB_trial_ec, save_path, s_name);
    num = num + 1;
    
end

%% PBA Data Extraction

% Get list of all available subject files
pba_files = get_files(pba_data_path, 'set');

% Set up for save name
num = 2000;

% Loop through all subjects, extracting and saving out data
for subj = pba_files
    
    s_name = [num2str(num), '.mat'];
    extract_save_data(subj, pba_data_path, pba_rest_ec, pba_trial_ec, save_path, s_name);
    num = num + 1;
    
end

