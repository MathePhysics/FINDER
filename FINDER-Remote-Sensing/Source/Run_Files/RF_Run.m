%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smooth_radar_data
% main script for smoothing radar data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%paths_hmm

methods.init = @RF_Initialize;%@init_smooth_radar_data;
methods.load_radar_data = @load_radar_data;
methods.smooth = @smooth_using_temporal_bayesian_filter_pcg_big_datastructure_pre;%@smooth_using_temporal_bayesian_filter;

parameters = methods.init();
radar_data = methods.load_radar_data(parameters);
output = methods.smooth(parameters,radar_data);

save(char(parameters.radar_output_folder+"filtered_radar_data_c_"+parameters.c+"_c2_"+parameters.c2+"_pcg_combined_all"+".mat"),'output','-v7.3')