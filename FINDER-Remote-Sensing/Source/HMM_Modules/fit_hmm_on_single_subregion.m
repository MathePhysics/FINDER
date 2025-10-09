%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fit_hmm_on_single_subregion
% Similar to fit_hmm: fits hidden markov model on a region. This function
% takes an array as an input instead of Output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subregion is a [sizeY,sizeX,number_of_days] array.

% TODO:
% - update code to use subregion instead of Output.
% - How to handle forest mask? How about we also include it as an input. 

function [subregion] = fit_hmm_on_single_subregion(parameters,methods,subregion,forestmask_subregion)


[height,width,n_days,~] = size(subregion);

Data_Queue = parallel.pool.DataQueue;
afterEach(Data_Queue, @print_progress);

temporary_combined_subregion_data = single(zeros(height,width,n_days));

rows_done = 0;

%% If there is a forest mask supplied, there is a check before fitting the HMM on any given pixel.

if ~isempty(parameters.forestmask_path)

    parfor (x = 1:width,1)    

                for y = 1:height
                if(forestmask_subregion(y,x) == 1)
                    pixel_labels = fit_hmm_on_single_subregion_pixel(y,x,parameters,subregion,methods);         
                    temporary_combined_subregion_data(y,x,:) = pixel_labels;

                end

                end

                rows_done = rows_done+1
  
                send(Data_Queue,x)
    end
        
 
  
    
else

    %% If there is no forest mask path supplied then there is no check; all pixels will be fit.
    
     parfor (x = 1:width,1)    

                for y = 1:height             

                    pixel_labels = fit_hmm_on_single_subregion_pixel(y,x,parameters,subregion,methods);   
                    temporary_combined_subregion_data(y,x,:) = pixel_labels;

                end

                rows_done = rows_done+1

                send(Data_Queue,x)
     end
     
       
end

  subregion(:,:,:,4) = temporary_combined_subregion_data;   


  function print_progress(~)

        rows_done = rows_done+1;        

  end

end