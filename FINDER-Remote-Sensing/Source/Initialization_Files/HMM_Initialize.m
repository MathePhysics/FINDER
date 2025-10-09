%%% Memory Initialize HMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set datafile and HMM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters = HMM_Initialize()
 % inputs: none
 % outputs: parameter data structure (data files/paths, transition
 % probabilities)

%1 for visual only,2 for radar only, 3 for hybrid
parameters.flag=3;

parameters.forestmask_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/forestmask.mat';
%parameters.optical_days='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/sentinel2_evi_days.mat';
parameters.processed_anomaly_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/KL/AnomalySeq_s2_full_maxtime_71_numEigen_29.mat'%'/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/Caleb_AnomalySeq_sentinel2_evi_small_area_1_full_maxtime_71_numEigen_29.mat';	
parameters.radar_data_path = '/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/RF/filtered_radar_data_c_1000_c2_10_pcg_combined_all.mat'%'/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/s1_combined_edited_smoothed_small_area_1_2022.mat';%"/projectnb/uqproj/Fusion/Caleb/Russell's Local/data/Radar/filtered_radar_data_c_1000_c2_10.mat";
parameters.radar_days_path ='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/s1_days.mat'%'/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/HMM/s1_combined_edited_days_366_2022.mat';%'/projectnb/uqproj/Fusion/data/Small_Area1/Sentinel1_days.mat';

parameters.HMM_output_file='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Output/HMM/la_hybrid.mat';
% -- radar parameters -----------------------------------------------------
parameters.radar_index = 1; 
% -------------------------------------------------------------------------

%%% HMM Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thresholds --------------------------------------------------------------
parameters.anomaly_threshold = 1.2; % more than this means likely anomaly
parameters.radar_threshold = -5.5; % less than this means likely bareground
%--------------------------------------------------------------------------
parameters.maxtime=71;
parameters.n_frames_to_confirm_class = 10;
%--------------------------------------------------------------------------
%% Subregion options
parameters.subregion_height = 1000;%461;
parameters.subregion_width  = 1000;%611;

%--------------------------------------------------------------------------
%% Additional Options 
parameters.post_flag=1;
parameters.post_count=20;

parameters.zero_out_all_visual_data = false;
parameters.kept_visual_days = []%[3774,3779,3784,3789];
parameters.radar_days_skipped=0;

%% Parallel Pool Options

parameters.number_of_workers = 28;



% probabilities -----------------------------------------------------------

if(parameters.flag==1.)
    parameters.transition_matrix = [0.950,0.024375,0.000625, 0.025;
                     0.950,0.024375,0.000625, 0.025 ;
                     0.025, 0.000625,0.024375,0.95;
                     0.025, 0.000625,0.024375,0.95];  
    
    parameters.anomaly_emission_probabilities = [0.005, 0.99,0.99,0.99];
    parameters.anomaly_emission_probabilities = [parameters.anomaly_emission_probabilities; 1-parameters.anomaly_emission_probabilities];
    
    parameters.initial_probabilities = [1/3,1/6,1/6,1/3]; 
end

if(parameters.flag==2.)

    % states:
    % s1: forest
    % s2: bareground
    parameters.transition_matrix = [0.95+0.024375,0.000625+0.025;
                                    0.025+0.000625,0.95+0.024375];
                
    parameters.radar_emission_probabilities = [0.35,0.8];
    parameters.radar_emission_probabilities = [parameters.radar_emission_probabilities; 1- parameters.radar_emission_probabilities];
    

    parameters.initial_probabilities = [1/2,1/2];   
end

if(parameters.flag==3.)
    % states:
    % s1: forest+no cloud
    % s2: forest+cloud
    % s3: bareground+cloud
    % s4: bareground+no cloud
    
    parameters.transition_matrix = [0.950,0.024375,0.000625, 0.025;
                     0.950,0.024375,0.000625, 0.025 ;
                     0.025, 0.000625,0.024375,0.95;
                     0.025, 0.000625,0.024375,0.95];             
                
    parameters.radar_emission_probabilities = [0.35,0.35,0.8,0.8];
    parameters.radar_emission_probabilities = [parameters.radar_emission_probabilities; 1- parameters.radar_emission_probabilities];
    
    parameters.anomaly_emission_probabilities = [0.005, 0.99,0.99,0.99];
    parameters.anomaly_emission_probabilities = [parameters.anomaly_emission_probabilities; 1-parameters.anomaly_emission_probabilities];
    parameters.combined_emission_probabilities = [parameters.anomaly_emission_probabilities(1,:).*parameters.radar_emission_probabilities(1,:); % both anomalous
                                 parameters.anomaly_emission_probabilities(2,:).*parameters.radar_emission_probabilities(1,:); % only radio anomalous
                                 parameters.anomaly_emission_probabilities(1,:).*parameters.radar_emission_probabilities(2,:); % only visual anomalous
                                 parameters.anomaly_emission_probabilities(2,:).*parameters.radar_emission_probabilities(2,:)]; % neither anomalous
    parameters.initial_probabilities = [1/3,1/6,1/6,1/3]; 
end

end
