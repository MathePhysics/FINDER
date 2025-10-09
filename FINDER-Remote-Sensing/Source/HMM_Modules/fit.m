function [bitmap_and_datemap] = fit(height_indices,width_indices,visual,radar,forestmask,parameters,methods)

[subregion_radar,subregion_visual,forestmask_subregion] = subregion_repackager(height_indices,width_indices,visual,radar,forestmask,parameters);
[subregion,parameters]=preprocess(parameters,subregion_visual,subregion_radar,length(height_indices),length(width_indices));

processed_data = fit_hmm_on_single_subregion(parameters,methods,subregion.combined_data,forestmask_subregion);
processed_data_struct=struct('combined_data',processed_data);
       
[accumulated_data,~] = methods.accumulate_classes(parameters,processed_data_struct);

bitmap_and_datemap=make_bitmap_and_datemap(accumulated_data,parameters.flag);

end