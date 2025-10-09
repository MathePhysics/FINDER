%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load_radar_data
% load radar data according to parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [radar_data] = load_radar_data(parameters)

tic
fprintf("Loading radar data...\n")
radar_data = load(parameters.radar_data_path);
fprintf("Finished loading radar data.\n")
toc

end