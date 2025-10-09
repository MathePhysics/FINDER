restoredefaultpath % removes everything from path 

% get the filepath of *this file* and change the directory to it
filepath_of_this_file = mfilename('fullpath');
filepath_of_this_file = filepath_of_this_file(1:(end-9)); %(-8 is to remove /paths.m)

%cd(filepath_of_this_file)

% adds only this pipeline folder and subfolders.
addpath(genpath("/projectnb/multifusion/Projects/Code/Code/Organized_Code")); 
