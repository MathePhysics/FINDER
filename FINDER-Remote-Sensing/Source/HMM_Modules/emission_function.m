%%% emission_function_hmm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function for computing probability of emissions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function emission_probability = emission_function(anomaly_value,radar_value,state,parameters,x,y)

% Return the probability of a given emission conditional on what is
% observed
% anomaly_value and radar_value are the observations
% state is the hidden underlying state index from 1-4: Forest,
% Forest+Cloud, Bareground_Cloud, Bareground


% Emission probabilities -----------------------------------------
% -----------------------------------------------------------------

if ~isnan(anomaly_value) && ~isnan(radar_value)  
    

    emission_probability = parameters.combined_emission_probabilities([true, false,true,false] == anomaly_value & [true, true, false, false] == radar_value,state);
elseif ~isnan(anomaly_value) && isnan(radar_value)

    emission_probability = parameters.anomaly_emission_probabilities(2-anomaly_value,state);

elseif isnan(anomaly_value) && ~isnan(radar_value)
    emission_probability = parameters.radar_emission_probabilities(2-radar_value,state);

elseif isnan(anomaly_value) && isnan(radar_value)

    emission_probability = 1;
    %fprintf("Warning: there should be no index missing both anomaly and radar")

end

end