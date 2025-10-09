function array = CompPredictAUC(Datas, parameters, methods)
    
% Predict SVM

switch parameters.multilevel.svmonly
    case 1
        SVMField = 'multilevel';
    case 0
        SVMField = 'multilevel';
    case 2
        SVMField = 'multilevel';
end


[class_Test_A,y_Test_A] = methods.all.SVMpredict(parameters.(SVMField).SVMModel, Datas.X_Test_A);
[class_Test_B,y_Test_B] = methods.all.SVMpredict(parameters.(SVMField).SVMModel, Datas.X_Test_B);

NA = size(Datas.X_Test_A, 1); NB = size(Datas.X_Test_B,1);

scores_A = nan(NA,1); scores_B = nan(NB,1);

scores_A(class_Test_A == 0) = min( y_Test_A(class_Test_A == 0,:), [], 2 );
scores_A(class_Test_A == 1) = max( y_Test_A(class_Test_A == 1,:), [], 2 );
scores_B(class_Test_B == 0) = min( y_Test_B(class_Test_B == 0,:), [], 2 );
scores_B(class_Test_B == 1) = max( y_Test_B(class_Test_B == 1,:), [], 2 );


labels = [ones(NA, 1); zeros(NB, 1)];
scores = [scores_A(:) ; scores_B(:)];



%[Xroc,Yroc] = perfcurve(labels, scores, 1);
switch parameters.data.validationType
    case 'Synthetic', nY = 2*parameters.synthetic.NTest;
    case 'Cross', nY = parameters.cross.NTestA + parameters.cross.NTestB;
    case 'Kfold', nY = 2*parameters.Kfold;
end
numpad = nY - length(scores);

if numpad > 0
    padnan = @(x) padarray(x(:), [numpad, 0], nan, 'post');
    scores= padnan(scores);
    labels = padnan(labels);
end

array = [labels(:),scores(:)] ;
array = reshape(array, [1,1,size(array)]);






end
