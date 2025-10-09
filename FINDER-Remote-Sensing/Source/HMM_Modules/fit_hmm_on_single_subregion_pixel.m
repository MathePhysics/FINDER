%%% fits a hidden markov model at a single pixel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pixel_labels = fit_hmm_on_single_subregion_pixel(y,x,parameters,subregion,methods)
% pixel labels is a vector with length equal to the number of total
% observations.
            combined_data = subregion;
            P = parameters.transition_matrix;

            anomaly_sequence = abs(squeeze(combined_data(y,x,:,1)));  
            radar_sequence = squeeze(combined_data(y,x,:,2)); 
            
            if(parameters.flag==3)
            both_nan_indices = isnan(anomaly_sequence) & isnan(radar_sequence); 
            anomaly_sequence = anomaly_sequence(~both_nan_indices);
            radar_sequence = radar_sequence(~both_nan_indices);
            end

            if(parameters.flag==1)
            nan_indices=isnan(anomaly_sequence);
            anomaly_sequence=anomaly_sequence(~nan_indices);
            end

            n=0;
            if(parameters.flag~=2)
            n = length(anomaly_sequence);
            else
            n=parameters.number_of_observation_days;
            end  
            
            if(n~=0)
            % Viterbi Algorithm ------------------------------------------------
            if(parameters.flag~=2)
            v = NaN(4,n);
            else
            v=NaN(2,n);
            end
            %Initialization step:
            if(parameters.flag==3)
            v(:,1) = parameters.initial_probabilities .* methods.emission_function(anomaly_sequence(1),radar_sequence(1),1:4,parameters,x,y);
            end
            if(parameters.flag==1)
            v(:,1) = parameters.initial_probabilities .* parameters.anomaly_emission_probabilities(2-anomaly_sequence(1),1:4);%emission_function_hmm_visual_only(anomaly_sequence(1),1:4,parameters);
            end
            if(parameters.flag==2)
            v(:,1) = parameters.initial_probabilities .* parameters.radar_emission_probabilities(2-radar_sequence(1),1:2); %emission_function_hmm_radar_only(radar_sequence(1),1:2,parameters);
            end
            %Induction step:
            for t = 2:n    
               for k = 1:4
                   if(parameters.flag==3)
                   v(k,t) = max(v(:,t-1) .* P(:,k) .* methods.emission_function(anomaly_sequence(t),radar_sequence(t),k,parameters,x,y));
                   end
                   if(parameters.flag==1)
                   %disp(anomaly_sequence(t));
                   v(k,t) = max(v(:,t-1) .* P(:,k) .* parameters.anomaly_emission_probabilities(2-anomaly_sequence(t),k));%emission_function_hmm_visual_only(anomaly_sequence(t),k,parameters));
                   %disp('ccc');
                   end 
                   if(parameters.flag==2 && k<=2)
                   v(k,t) = max(v(:,t-1) .* P(:,k) .* parameters.radar_emission_probabilities(2-radar_sequence(t),k));
                   end   

                       
               end

            end
            %Extract Sequence:
            %disp('ddd');
            pixel_labels = nan(n,1); % this will label 1 = normal, 2 = cloud, 3 = anomaly
            
            [~,pixel_labels(n)] = max(v(:,n));

            for t = (n-1):(-1):1
                [~,pixel_labels(t)] = max(v(:,t) .* P(:,pixel_labels(t+1))); 
            end
            if(parameters.flag==3)
            pixel_labels(pixel_labels == 3) = 2; % cloud
            pixel_labels(pixel_labels == 4) = 3; % deforestation     
            
            vec_holder=nan(length(both_nan_indices),1);
            vec_holder(~both_nan_indices)=pixel_labels;
            pixel_labels=vec_holder;
            end

            if(parameters.flag==1)
            vec_holder=nan(length(nan_indices),1);
            vec_holder(~nan_indices)=pixel_labels;
            pixel_labels=vec_holder;
            end

            else
            pixel_labels=nan(length(nan_indices),1);
            end
            
end

            

