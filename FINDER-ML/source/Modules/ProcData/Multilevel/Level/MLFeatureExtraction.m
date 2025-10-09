function [X,y] = MLFeatureExtraction(parameters, l, value, class)

%C = parameters.Training.A.C;
C = parameters.(value).(class).C(:,1:end-1);

m = size(C,2);
% nested choosed 0~l level
%lastColumn = C(:, end);
lastColumn = parameters.(value).(class).levelcoeff;

switch parameters.multilevel.nested
    case 0 
        X = C(lastColumn == l, :);
    case 1
        X = C(lastColumn <= l & lastColumn >= 0, :); 
    case 2
        X = C(lastColumn >= l, :);
end

%lastColumn = X(:, end);
%X = X(lastColumn >=0, :);

%X(:,end)=[];  %delete the level column

%y = ones(1, m-1);
switch class
    case 'A', y = ones(1,m);
    case 'B', y = zeros(1,m);
end

X = X';
y = y';

%Check:
% X2 = parameters.Training.A.C(parameters.Training.A.levelcoeff >= l,1:end-1);
% disp(all(all(X2' == X)));
% stem(lastColumn >= l)
% hold on

end