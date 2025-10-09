function Datas = MinimizeExpectedAngle(Datas, parameters, methods)

if ~parameters.transform.ComputeTransform, return, end

%RA = min(parameters.data.numofgene, parameters.data.A);
[H0, ~, ~] = svd([Datas.A.CovTraining Datas.B.CovTraining], 'econ');
H00 = H0; %Save copy for comparison later
%H0 = speye(parameters.data.numofgene, RA);
%H0 = speye(parameters.data.numofgene, RA);
%H0 = randn(parameters.data.numofgene, RA);
itercounter = 0;
nincrease = 0;
BestObj = Inf;
obj0 = 1;
normgrad = Inf;
objList = [];
BestObjList = [];
stepdiff = Inf;
keepIterating = true;
maxiter = 250; 
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
    Sjk = HSIP(H0, H0, AAll, BAll);
    Sjj = HSN(H0, AAll);
    Skk = HSN(H0, BAll);
    Ljk = (Sjk.^2) ./ (Sjj(:) .* Skk(:)');
    objn = mean(Ljk, 'all' );

    
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
    fprintf('Objective Value: %0.2f. ', cos2deg(objn));
    fprintf('Best Objective: %0.2f. ', cos2deg(BestObj));
    fprintf('Relative Change: %0.2e.', stepdiff);
    fprintf('Gradient Norm: %0.2e \n', normgrad)
    end

    %if objn > obj0, nincrease = nincrease + 1; end

%% Decide to continue
    x = [objn < eps, % Objective is already close to zero
         stepdiff < eps, % Relative decrease in objective function is close to zero
         itercounter > maxiter, % Exceed the maximum number of iterations 
         normgrad < eps, %Gradient is close to zero
         nincrease >= maxincrease];  %Objective function increased too many times
    
    if any(x)
        break
    end

    obj0 = objn;

%% Estimate Gradient
    Sjk = HSIP(H0, H0, ATrain, BTrain);
    Sjj = HSN(H0, ATrain);
    Skk = HSN(H0, BTrain);
    Ljk = (Sjk.^2) ./ (Sjj(:) .* Skk(:)');
    Lj = sum(Ljk,2); Lk = sum(Ljk,1);
    Lj = Lj(:)'; Lk = Lk(:)';


    B(:,:,1) =  1 * ATrain * (Ljk./Sjk) * (BTrain' * H0) ; 
    B(:,:,2) =  1 * BTrain * (Ljk./Sjk)' * (ATrain' * H0) ;
    B(:,:,3) =  -1 * (ATrain .* (Lj./Sjj)) * (ATrain' * H0) ; 
    B(:,:,4) =  -1 * (BTrain .* (Lk./Skk)) * (BTrain' * H0) ; 

    grad = sum(B,3) /(NATrain * NBTrain) ;
  
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
    t = 1/(sqrt(itercounter) * normgrad);
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
fprintf('Difference in Final vs Initial H: %0.2f \n\n', norm(BestH - H00, 'fro'));

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