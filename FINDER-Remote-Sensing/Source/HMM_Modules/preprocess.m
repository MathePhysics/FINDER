%%% preprocess_hmm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -load & transform data based on hmm parameters
%%% -categorize data based on thresholds & categorize available observations
%%% (radar, visual, or both)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Output,parameters] = preprocess(parameters,visual,radar,subregion_height,subregion_width)



%% Load Data ---------------------------------------------------------------
% -Imports both radar and visual data
% -If no path for radar data is supplied, radar will be set to NaN
% -variable maxtime stores last training frame
% -the anomaly data loaded is assumed to contain both the training and
% testing slides, hence the maxtime argument is used. This is probably not
% the case.

fprintf("Loading Data...")

Output.anomaly_data=visual;
Output.raw_radar_data = radar;

parameters.observation_data.days_with_visual_observation=[];
parameters.observation_data.days_with_radar_observation=[];

if(parameters.flag==1)    
    Output.raw_radar_data = [];
    processed_radar_data = NaN;
end    

if(parameters.flag==2)    
    Output.anomaly_data = [];
    processed_anomaly_data = NaN;
end  

if(parameters.flag~=2)
    % Process Anomaly Data. Entries that are NaN should stay NaN.
    maxtime = parameters.maxtime;%Output.anomaly_data.parameters.data.maxtime;

    missing_indices = isnan(Output.anomaly_data.collectresults);
    processed_anomaly_data = 1.*(abs(Output.anomaly_data.collectresults) > parameters.anomaly_threshold);

    
    processed_anomaly_data(missing_indices) = NaN;
    

    processed_anomaly_data = abs(processed_anomaly_data(:,:,(maxtime+1):end)); % Removes training frames
    try
        parameters.observation_data.days_with_visual_observation = Output.anomaly_data.Output.info_days((maxtime+1):end);  
    catch
        try
            parameters.observation_data.days_with_visual_observation = Output.anomaly_data.parameters.data.dias(72:end); %Output.anomaly_data.Output.info_days((maxtime+1):end);  
        catch
            parameters.observation_data.days_with_visual_observation = Output.anomaly_data.parameters.data.dias(maxtime+1:end); 
        end
    end

end


if(parameters.flag~=1)
    skipped=parameters.radar_days_skipped;
    processed_radar_data = Output.raw_radar_data(:,:,skipped+1:end,parameters.radar_index);
    missing_indices = isnan(processed_radar_data);
    
    processed_radar_data = 1.*(processed_radar_data < parameters.radar_threshold);
    processed_radar_data(missing_indices) = NaN;
    skipped=parameters.radar_days_skipped;
    parameters.observation_data.days_with_radar_observation = parameters.radar_days(skipped+1:end);%Output.raw_radar_data(1,1,skipped+1:end,2);
    

end

parameters.observation_data.days_with_any_observation = union(parameters.observation_data.days_with_visual_observation,parameters.observation_data.days_with_radar_observation);
parameters.number_of_observation_days = length(parameters.observation_data.days_with_any_observation);

%% POSSIBLE UPDATE: [[Import anomaly data first so that it isn't loaded twice when there is no radar data.]]

% If there is no radar data, the radar data is set to NaN.


% raw_radar_data is a 4D array: 
% -X x Y x Day x Datatype. Datatype 1&2 are date information, 3&4 radar.
%;
%end
% raw_anomaly_data is a structure with Output/Parameters
% - raw_anomaly_data.Ouptupt.collectresults is the anomaly values
% - raw_anomaly_data.parameters.data.patchsizeI and J are patch locations


%maxtime = Output.anomaly_data.parameters.data.maxtime; % last training frame. 

fprintf("DONE.  \n")



%% process anomaly and radar data using thresholds ------------------------

fprintf("Processing Data...")






% generate combined dataset

Output.combined_data = NaN(subregion_height,subregion_width,parameters.number_of_observation_days,5); 

Output.combined_data(:,:,ismember(parameters.observation_data.days_with_any_observation,parameters.observation_data.days_with_visual_observation),1) = processed_anomaly_data;
Output.combined_data(:,:,ismember(parameters.observation_data.days_with_any_observation,parameters.observation_data.days_with_radar_observation),2) = processed_radar_data;


% include day information in combined data
for ii = 1:subregion_height
    for jj = 1:subregion_width
        Output.combined_data(ii,jj,:,3) = parameters.observation_data.days_with_any_observation;
    end
end                        
              

if ~isempty(parameters.kept_visual_days) || parameters.zero_out_all_visual_data
    
    %combined_data_days = parameters.observation_data.days_with_any_observation;        
   % Output.combined_data(:,:,~ismember(combined_data_days,parameters.kept_visual_days),1) = NaN;
    %parameters.observation_data.days_with_visual_observation = parameters.kept_visual_days;    

end






%fprintf("Determining pixelwise observation availability...")

% 1 = both, 2 = anomaly only, 3 = radar only, 4 = neither

%Output.combined_data(:,:,:,4) = ismissing(Output.combined_data(:,:,:,1))+is.missing(Output.combined_data(:,:,:,2))

%fprintf("DONE.")


fprintf("DONE.  \n")
end