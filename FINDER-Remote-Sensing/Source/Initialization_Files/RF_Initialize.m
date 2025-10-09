%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% init_smooth_radar_data
% Generate parameters for smooth_radar_data script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parameters] = RF_Initialize()

parameters.radar_data_path = '/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/RF/s1.mat';
parameters.radar_output_folder='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/RF/';

parameters.L = 1;
parameters.c = 1000; % spatial smoothing penalty
parameters.c2 = 10; % temporal smoothing penalty
parameters.pcg.maxit = 200; % for pcg
parameters.pcg.tol = 1e-6;

end