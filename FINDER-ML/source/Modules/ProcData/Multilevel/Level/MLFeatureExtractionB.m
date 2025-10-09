function [X, y] = MLFeatureExtractionB(parameters, l)

%C = parameters.Training.B.C;
C = parameters.Training.B.C(:,1:end-1);

m = size(C,2);

% nested choosed 0~l level
%lastColumn = C(:, end);
%X = C(lastColumn <= l, :); 
lastColumn = parameters.Training.B.levelcoeff;

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

%y = zeros(1, m-1);
y = zeros(1,m);
X = X';
y = y';

%Check
% X2 = parameters.Training.B.C(parameters.Training.B.levelcoeff >= l,1:end-1);
% disp(all(all(X2' == X)));
% %stem(X2(:,1))


end