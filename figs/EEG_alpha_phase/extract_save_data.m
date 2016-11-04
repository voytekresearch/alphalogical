function extract_save_data(subj, path, rest_ec, trial_ec, s_path, s_name)
%%
%

%%
data = pop_loadset('filename', subj, 'filepath', path);

rest_data = pop_epoch(data, rest_ec, [4, 124]);
oz_rest_data = rest_data.data(30, :);

trial_data = pop_epoch(data, trial_ec, [0, 120]);
oz_trial_data = trial_data.data(30, :, 1);

ev_types = trial_data.epoch(1).eventtype;
ev_times = trial_data.epoch(1).eventlatency;

%%
save([s_path, s_name], 'oz_rest_data', 'oz_trial_data', 'ev_types', 'ev_times')