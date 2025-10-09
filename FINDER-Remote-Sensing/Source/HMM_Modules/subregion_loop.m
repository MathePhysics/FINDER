function [bitmap_and_datemap] = subregion_loop(parameters,methods)


visual=NaN;
radar=NaN;
forestmask=NaN;

if ~isempty(parameters.forestmask_path)
    tic;
    fprintf("Loading Forest Mask...")
    forestmask = load(parameters.forestmask_path);
    fprintf("DONE \n");
    toc;

else
   
    fprintf("No Forest Mask path supplied: all pixels will be classified \n")
    
end

subregion_height = parameters.subregion_height;
subregion_width  = parameters.subregion_width;           


if(parameters.flag~=2)
    visual=matfile(parameters.processed_anomaly_path);

    try
        height=size(visual,'collectresults',1);
        width=size(visual,'collectresults',2);
    catch
        visual=load(parameters.processed_anomaly_path);
        vis_Output=visual.Output;
        vis_collectresults=visual.Output.collectresults;
        vis_parameters=visual.parameters;
        visual=struct('Output',vis_Output,'parameters',vis_parameters,'collectresults',vis_collectresults);
        height=size(visual.collectresults,1);
        width=size(visual.collectresults,2);  
    end
end

if(parameters.flag~=1)
        radar=matfile(parameters.radar_data_path);
        radar_days=load(parameters.radar_days_path);
        parameters.radar_days=radar_days.days;
    
        height=size(radar,'output',1);
        width=size(radar,'output',2);     
end

height_bounds = [1:(subregion_height):(height-subregion_height +1); subregion_height:(subregion_height):height];
width_bounds  = [1:(subregion_width):(width-subregion_width +1); subregion_width :(subregion_width ):width];

[~,n_vertical_subregions]   = size(height_bounds);
[~,n_horizontal_subregions] = size(width_bounds);
datemap=NaN(height,width);
bitmap=NaN(height,width);
% =========================== Run HMM on each subregion  ==================

fprintf("Fitting hmm on "+height+"x"+width+" grid of subregions...\n")

for v = 1:n_vertical_subregions
   for h = 1:n_horizontal_subregions    
 
       tic;
       height_indices = height_bounds(1,v):height_bounds(2,v);
       width_indices = width_bounds(1,h):width_bounds(2,h);
        
       bitmap_and_datemap = fit(height_indices,width_indices,visual,radar,forestmask,parameters,methods);
       datemap(height_indices,width_indices)=bitmap_and_datemap.datemap;
       bitmap(height_indices,width_indices)=bitmap_and_datemap.bitmap;
       
       fprintf("Model fitted on subregion ["+v+","+h+"]\n");
       toc;
   end    
      if(mod(width,subregion_width)~=0)
           tic;
           height_indices=height_bounds(1,v):height_bounds(2,v);
           width_indices=width_bounds(2,n_horizontal_subregions)+1:width;
    
           bitmap_and_datemap = fit(height_indices,width_indices,visual,radar,forestmask,parameters,methods); 
           datemap(height_indices,width_indices)=bitmap_and_datemap.datemap;
           bitmap(height_indices,width_indices)=bitmap_and_datemap.bitmap;
           
           h2=h+1;
           fprintf("Model fitted on edge subregion ["+v+","+h2+"]\n");
           toc;               
     end
end


if(mod(height,subregion_height)~=0)
    for h = 1:n_horizontal_subregions
       tic;
       height_indices=height-subregion_height+1:height;
       width_indices=(h-1)*subregion_width+1:h*subregion_width;
       
       bitmap_and_datemap = fit(height_indices,width_indices,visual,radar,forestmask,parameters,methods); 
        
       datemap(height_indices,width_indices)=bitmap_and_datemap.datemap;
       bitmap(height_indices,width_indices)=bitmap_and_datemap.bitmap;
       
       v2=v+1;
       fprintf("Model fitted on edge subregion ["+v2+","+h+"]\n");
       toc;        
    end
end
if(mod(height,subregion_height)~=0 & mod(width,subregion_width)~=0)
    tic;
    height_indices=height-subregion_height+1:height;
    width_indices=width-mod(width,subregion_width)+1:width;

    bitmap_and_datemap = fit(height_indices,width_indices,visual,radar,forestmask,parameters,methods);       
    datemap(height_indices,width_indices)=bitmap_and_datemap.datemap;
    bitmap(height_indices,width_indices)=bitmap_and_datemap.bitmap;
     
    v2=v+1;
    h2=h+1;
    fprintf("Model fitted on corner subregion ["+v2+","+h2+"]\n");
    toc;          
end
if(parameters.post_flag==1)
    datemap=postprocess(datemap,parameters.post_count);
end  

bitmap_and_datemap=struct('bitmap',bitmap,'datemap',datemap);
end