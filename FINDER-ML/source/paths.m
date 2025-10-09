d = fileparts(pwd);
a = genpath(d);
path(path,a);

x = dir(pwd);
isML_Features = arrayfun(@(y) contains(y.name, 'ML_Features'), x);
ML_FeatureFiles = x(isML_Features);
arrayfun(@(y) delete(y.name), ML_FeatureFiles);
