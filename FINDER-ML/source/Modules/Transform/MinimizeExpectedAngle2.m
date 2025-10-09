function Datas = MinimizeExpectedAngle(Datas, parameters, methods)

if ~parameters.transform.ComputeTransform, return, end

RA = min(parameters.data.numofgene, parameters.data.A);
[H0, ~, ~] = svd(Datas.A.CovTraining, 'econ');
%H0 = speye(parameters.data.numofgene, RA);
%H0 = speye(parameters.data.numofgene, RA);
%H0 = randn(parameters.data.numofgene, RA);
itercounter = 0;
nincrease = 0;
BestObj = Inf;
obj0 = 1;
objList = [];
BestObjList = [];
stepdiff = Inf;
keepIterating = true;
maxiter = 250; 
normgrad = Inf;
maxincrease = 5;

HSIP = @(X,H, A,B) (A'*X) * (H'*B) ;
HSIP2 = @(X,H, S, A, B) (X'*A) * S * (B'*H); 
HSN = @(H,A) sum( (H'*A).^2, 1);
cos2deg = @(x) acos(sqrt(x)) * 180 /pi;


weights = [1 / (parameters.data.A * parameters.data.B)
          1 / (parameters.data.A * parameters.data.B)
          -1 / (parameters.data.A )
           -1 / (parameters.data.B) ];



%% Leave out sample for testing;
ATrain = Datas.A.CovTraining; 
BTrain = Datas.B.CovTraining; 
%BTrain = Datas.B.CovTraining - mean(Datas.B.Training,2); 

if ~strcmp(parameters.data.validationType, 'Cross')
        ATest = ATrain(:,1:parameters.Kfold);
        BTest = BTrain(:,1:parameters.Kfold);
        ATrain = ATrain(:,parameters.Kfold+1:end);
        BTrain = BTrain(:,parameters.Kfold+1:end);
else
        ATest = ATrain(:,1:parameters.cross.NTestA);
        BTest = BTrain(:,1:parameters.cross.NTestB);
        ATrain = ATrain(:,parameters.cross.NTestA+1:end);
        BTrain = BTrain(:,parameters.cross.NTestB+1:end);
end

AAll = [ATrain ATest];
BAll = [BTrain BTest];
NATrain = size(ATrain,2);
NBTrain = size(BTrain,2);


while keepIterating
%% Increment iteration
itercounter = itercounter + 1;
 

%% Estimate objective
    SAB_jk = HSIP(H0, H0, AAll, BAll);
    SAA_jj = HSN(H0, AAll);
    SBB_kk = HSN(H0, BAll);
    SAA_jk = HSIP(H0, H0, AAll, AAll);
    SBB_jk = HSIP(H0, H0, BAll, BAll);

    LAB = (SAB_jk.^2) ./ (SAA_jj(:) .* SBB_kk(:)');
    LAA = (SAA_jk.^2) ./ (SAA_jj(:) .* SAA_jj(:)');
    LBB = (SBB_jk.^2) ./ (SBB_kk(:) .* SBB_kk(:)');

    avgAngle = cellfun( @(x) mean(x, 'all'), {LAB, -LAA, -LBB});
    objn = 2 + sum(avgAngle);

    
    stepdiff = abs((objn - obj0)/obj0);

    if itercounter == 1, InitObj = objn; end

    if objn < BestObj
        BestObj = objn;
        BestH = H0;
        %nincrease =0;
    elseif objn >= BestObj
        %nincrease = nincrease + 1;
    end 

    if objn < obj0, nincrease = 0; elseif objn > obj0, nincrease = nincrease + 1; end

    
    if mod(itercounter, floor(maxiter/10)) == 0
    fprintf('Iteration: %d \n', itercounter)
    fprintf('Objective Value: %0.2e. ', objn);
    fprintf('Best Objective: %0.2e', BestObj);
    fprintf('Relative Change: %0.2e \n', stepdiff);
    end

    %if objn > obj0, nincrease = nincrease + 1; end

%% Decide to continue
    x = [objn < -2, % Objective is already close to zero
         stepdiff < eps, % Relative decrease in objective function is close to zero
         itercounter > maxiter, % Exceed the maximum number of iterations 
         normgrad < eps, %Gradient is close to zero
         nincrease >= maxincrease];  %Objective function increased too many times
    
    if any(x)
        break
    end

    obj0 = objn;

%% Estimate Gradient

    %AB crossover
    SAB_jk = HSIP(H0, H0, ATrain, BTrain);
    SAA_jj = HSN(H0, ATrain);
    SBB_kk = HSN(H0, BTrain);
    LAB = (SAB_jk.^2) ./ (SAA_jj(:) .* SBB_kk(:)');
    LABj = sum(LAB,2); LABk = sum(LAB,1);
    LABj = LABj(:)'; LABk = LABk(:)';


    B(:,:,1) =  1 * ATrain * (LAB./SAB_jk) * (BTrain' * H0) ; 
    B(:,:,2) =  1 * BTrain * (LAB./SAB_jk)' * (ATrain' * H0) ;
    B(:,:,3) =  -1 * (ATrain .* (LABj./SAA_jj)) * (ATrain' * H0) ; 
    B(:,:,4) =  -1 * (BTrain .* (LABk./SBB_kk)) * (BTrain' * H0) ; 

    %A crossover
    SAA_jk = HSIP(H0, H0, ATrain, ATrain);
    LAA = (SAA_jk.^2) ./ (SAA_jj(:) .* SAA_jj(:)');
    LAj = sum(LAA,2); LAk = sum(LAA,1);
    LAj = LAj(:)'; LAk = LAk(:)';

    B(:,:,5) =  1 * ATrain * (LAA./SAA_jk) * (ATrain' * H0) ; 
    B(:,:,6) =  1 * ATrain * (LAA./SAA_jk)' * (ATrain' * H0) ;
    B(:,:,7) =  -1 * (ATrain .* (LAj./SAA_jj)) * (ATrain' * H0) ; 
    B(:,:,8) =  -1 * (ATrain .* (LAk./SAA_jj)) * (ATrain' * H0) ; 

    %B Crossover
    SBB_jk = HSIP(H0, H0, BTrain, BTrain);
    LBB = (SBB_jk.^2) ./ (SBB_kk(:) .* SBB_kk(:)');
    LBj = sum(LBB,2); LBk = sum(LBB,1);
    LBj = LBj(:)'; LBk = LBk(:)';


    B(:,:,9) =  1 * BTrain * (LBB./SBB_jk) * (BTrain' * H0) ; 
    B(:,:,10) =  1 * BTrain * (LBB./SBB_jk)' * (BTrain' * H0) ;
    B(:,:,11) =  -1 * (BTrain .* (LBj./SBB_kk)) * (BTrain' * H0) ; 
    B(:,:,12) =  -1 * (BTrain .* (LBk./SBB_kk)) * (BTrain' * H0) ; 

    
    grad = sum(B(:,:,1:4),3) /(NATrain * NBTrain) -...
           sum(B(:,:,5:8),3) / (NATrain^2) - ...
           sum(B(:,:,9:12),3) / (NBTrain^2);
 %% Estimate Hessian 
    
    
    % C{1} = HSIP2(grad, grad, 1./S{1}, CovATrain, CovBTrain) + ...
    %        HSIP2(grad, grad, 1./S{2}, CovBTrain, CovATrain)';
    % C{2} = diag(HSIP2(grad, grad, diag(1./S{3}), CovATrain, CovATrain).^2);
    % C{3} = diag(HSIP2(grad, grad, diag(1./S{4}), CovBTrain, CovBTrain).^2);
    % 
    % Hess1  = cellfun(@(C) mean((C), 'all'), C);
    % Hess1 = [2 -1 -1]*Hess1(:);
    % 
    % D{1} = HSIP2(grad, H0, 1./S{1}, CovATrain, CovBTrain).^2 + ...
    %        HSIP2(grad, H0, 1./S{2}, CovBTrain, CovATrain).^2 ;
    % D{2} = diag(HSIP2(grad, H0, diag(1./S{3}), CovATrain, CovATrain).^2);
    % D{3} = diag(HSIP2(grad, H0, diag(1./S{4}), CovBTrain, CovBTrain).^2);
    % 
    % Hess2  = cellfun(@(C) mean((C), 'all'), C);
    % Hess2 = [-1 2 2]*Hess2(:);


 %% Compute Learning Rate
    % b = norm(grad,'fro')^2;
    % a = Hess1 + Hess2;
    % t = abs(b/a); 
    %t = norm(grad, 'fro') / itercounter;

%% Compute Learning Rate
    normgrad = norm(grad, 'fro');
    t = 1/(itercounter * normgrad );
     %t = 1 / norm(grad, 'fro');

 %% Compute step 
    H1 = H0 - t*grad;
    H1 = H1 / norm(H1, 'fro'); 
    %[U,S,~] = svd(H1, 'econ', 'vector');
    %H1 = U.*S';
    H0 = H1;
    %stepdiff = norm(H1 - H0, 'fro');
    
end



%fprintf(['\n ========== \n ' ...
%    'Initial Objective: %0.2f \n' ...
%    'Best Objective: %0.2f ' ...
%    '\n =========== \n'], cos2deg(InitObj), cos2deg(BestObj))


[U,S,V] = svd(BestH, 'econ', 'vector');
fprintf('Condition Number; %0.3e \n', max(S) / min(S) );

initTag = ["Initial", "Final"];

%Datas.A.All = Datas.rawdata.AData - mean(Datas.rawdata.AData,2);
%Datas.B.All = Datas.rawdata.BData;
for C = 'AB', Datas.(C).All = [Datas.(C).Training, Datas.(C).Testing]; end

Sets = ["CovTraining", "Testing", "Machine", "All"];

for i = 1:length(initTag)
 fprintf('%s Angles: \n ================= \n', initTag(i))
 tableArgs = {};
for Set = Sets
   
     if initTag(i) == "Final", for C = 'AB'
           Datas.(C).(Set) = (U .* S') * (U' * Datas.(C).(Set));
     end, end



SAB_jk = Datas.A.(Set)' * Datas.B.(Set);
SAA_jk = Datas.A.(Set)' * Datas.A.(Set);
SBB_jk = Datas.B.(Set)' * Datas.B.(Set);
SAA_jj = sum(Datas.A.(Set).^2, 1);
SBB_kk = sum(Datas.B.(Set).^2, 1);

LAB = (SAB_jk.^2) ./ (SAA_jj(:) .* SBB_kk(:)');
LAA = (SAA_jk.^2) ./ (SAA_jj(:) .* SAA_jj(:)');
LBB = (SBB_jk.^2) ./ (SBB_kk(:) .* SBB_kk(:)');



%avgAngle = cellfun( @(x) real(cos2deg(mean(x, 'all'))), {LAB, LAA, LBB}); 
avgAngle = cellfun( @(x) real((mean(cos2deg(x), 'all'))), {LAB, LAA, LBB}); 
avgAngle = avgAngle(:);

tableArgs = [tableArgs, {avgAngle}];
% avg_angle = cellfun(...
%             @(x) sprintf('%0.2f', cos2deg(objn));
% 
% 
% end
% fprintf('====================== \n')   
end

T = table(tableArgs{:}, 'VariableNames', Sets, 'RowNames', ["AB", "AA", "BB"]);
disp(T)
end

for C = 'AB', Datas.(C) = rmfield(Datas.(C), "All"); end



% fprintf('Initial Angles: \n ================= \n')
% for Set = ["CovTraining", "Testing", "Machine"]
% ZB = Datas.B.(Set);
% %ZB = Datas.B.(Set) - mean(Datas.B.(Set));
% Sjk = Datas.A.(Set)' * ZB;
% Sjj = sum(Datas.A.(Set).^2, 1);
% Skk = sum(ZB.^2, 1);
% Ljk = (Sjk.^2) ./ (Sjj(:) .* Skk(:)');
% objn = mean(Ljk, 'all' );
% 
% fprintf('%s angle: %0.2f \n', Set, cos2deg(objn));
% end
% fprintf('====================== \n')
% 
% 
% 
% %% Display Initial Angles
% fprintf('Final Angles: \n ================= \n')
% for Set = ["CovTraining", "Testing", "Machine"]
%  for C = 'AB'
%        Datas.(C).(Set) = (U .* S') * (U' * Datas.(C).(Set));
% end
% ZB = Datas.B.(Set);
% %ZB = Datas.B.(Set) - mean(Datas.B.(Set));
% Sjk = Datas.A.(Set)' * ZB;
% Sjj = sum(Datas.A.(Set).^2, 1);
% Skk = sum(ZB.^2, 1);
% Ljk = (Sjk.^2) ./ (Sjj(:) .* Skk(:)');
% objn = mean(Ljk, 'all' );
% 
% fprintf('%s angle: %0.2f \n', Set, cos2deg(objn));
% end

end