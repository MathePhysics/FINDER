function [subregion_radar,visual_struct,forestmask_subregion] = subregion_repackager(height_indices,width_indices,visual,radar,forestmask,parameters)


        subregion_radar=NaN;
       visual_struct=NaN;

       if(parameters.flag~=2)
         subregion_collectresults = visual.collectresults(height_indices,width_indices,:);
         visual_struct_parameters=visual.parameters;
         visual_struct_output=visual.Output;
         visual_struct=struct('collectresults',subregion_collectresults,'parameters',visual_struct_parameters,'Output',visual_struct_output);
       end

       if(parameters.flag~=1)
         subregion_radar = radar.output(height_indices,width_indices,:,:);
         start=1;
         subregion_radar=subregion_radar(:,:,start:end,:);
       end

       try
            forestmask_subregion = forestmask.forestmask(height_indices,width_indices); 
       catch
           try
                forestmask_subregion = forestmask.forestmask_2018(height_indices,width_indices);
           catch
                forestmask_subregion = forestmask.forestmask_2020(height_indices,width_indices);
           end
       end
end