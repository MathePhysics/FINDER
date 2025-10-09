%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make_bitmap_and_datemap
% Converts Output data structure into a bitmap of *final* classifications
% and a datemap of the day in which each classification was made
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitmap_and_datemap] = make_bitmap_and_datemap(Output,flag)

[ydim,xdim,n_days,~] = size(Output.combined_data);

bitmap_and_datemap.bitmap = Output.combined_data(:,:,n_days,4);

bitmap_and_datemap.datemap = nan(ydim,xdim);

for jj = 1:ydim
   for kk = 1:xdim
      
       if (flag==3 && Output.combined_data(jj,kk,n_days,4) == 3)
            bitmap_and_datemap.datemap(jj,kk) = Output.combined_data(1,1,find(Output.combined_data(jj,kk,:,4) == 3,1),3);

       end 

       if(flag==2 && Output.combined_data(jj,kk,n_days,4) == 2)
            bitmap_and_datemap.datemap(jj,kk) = Output.combined_data(1,1,find(Output.combined_data(jj,kk,:,4) == 2,1),3);
       end

       if(flag==1 && Output.combined_data(jj,kk,n_days,4) == 2)
            bitmap_and_datemap.datemap(jj,kk) = Output.combined_data(1,1,find(Output.combined_data(jj,kk,:,4) == 2,1),3);
       end       
       
   end
end

end