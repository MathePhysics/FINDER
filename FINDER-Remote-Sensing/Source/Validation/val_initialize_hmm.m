function parameters = val_initialize_hmm(at, rt, ftc, flag)

parameters.flag=flag;
% -- paths ----------------------------------------------------------------

parameters.input_path='/projectnb/multifusion/Projects/Code/Organized_Code/Data/Input/Validation/';
parameters.processed_anomaly_path=[parameters.input_path,'val_optical_rand.mat'];
parameters.radar_data_path = [parameters.input_path,'val_radar.mat'];
parameters.radar_days_path = [parameters.input_path,'s1_days_T.mat'];
parameters.forestmask_path = [parameters.input_path,'val_forestmask.mat'];
% -------------------------------------------------------------------------

% -- radar parameters -----------------------------------------------------
parameters.radar_index = 1; % array index for radar data (3 or 4)
% -------------------------------------------------------------------------

%%% HMM Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thresholds --------------------------------------------------------------
parameters.anomaly_threshold = at; % more than this means likely anomaly
parameters.radar_threshold = rt; % less than this means likely bareground
%--------------------------------------------------------------------------
parameters.maxtime=71;
parameters.n_frames_to_confirm_class = ftc;
%--------------------------------------------------------------------------
%% Subregion options
parameters.subregion_height = 1;%461;
parameters.subregion_width  = 1000;%611;

%--------------------------------------------------------------------------
%% Additional Options 

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
