%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Accumulate Classes
% Function to transform hmm output using some kind of simple rule, e.g "If
% something is deforestation/bareground for at least 10 days in a row"

function [Output,parameters] = accumulate_classes(parameters,Output)

[x,y,d,t] = size(Output.combined_data);

new_combined_data = zeros(x,y,d,t);
days=Output.combined_data(1,1,:,3);

if(parameters.flag==3)
bareground_classes=[3];
reset_classes=[1];
anomaly_val=3;
a=[1,3,4];
end

if(parameters.flag==1)
bareground_classes=[3,4];
reset_classes=[1,2];
anomaly_val=2;
a=[2];
end

if(parameters.flag==2)
bareground_classes=[2];
reset_classes=[1];
anomaly_val=2;
a=[2];
end

for jj = 1:x
   for kk = 1:y

       classes = squeeze(Output.combined_data(jj,kk,:,4));
       longstanding_classes = classes;
       
       longstanding_bareground = 0;
       
       for ll = 1:d
          
           if ismember(classes(ll),bareground_classes) & (days(ll)>=1540)
               longstanding_bareground = longstanding_bareground+1;
               longstanding_classes(ll) = longstanding_bareground;
           end
           

           
           if ismember(classes(ll),reset_classes)
                longstanding_bareground = 0;
           end     
                      
       end   
       
       classes(cummax(longstanding_classes) > parameters.n_frames_to_confirm_class) = anomaly_val;
       classes(cummax(longstanding_classes) <= parameters.n_frames_to_confirm_class & ismember(classes,a)) = 1;
       
       new_combined_data(jj,kk,:,4) = classes;
              
   end    
end

Output.combined_data(:,:,:,4) = new_combined_data(:,:,:,4);
parameters = parameters;

end