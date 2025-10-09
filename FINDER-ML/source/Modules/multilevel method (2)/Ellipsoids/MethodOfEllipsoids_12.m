function parameters = MethodofEllipsoids_12(Datas, parameters, methods)

return

%% Prep data
for C = 'AB'
    parameters.data.(C) = size(Datas.rawdata.([C 'Data']), 2);
end
Datas = InScriptPrepData(Datas, parameters, methods);


%Get a list of truncation parameters for class A
Truncations = GetTruncations(parameters);

%Initialize Array of wrongly misplaced points
Record = Inf;
P = parameters.data.numofgene;
M = length(Truncations);
N = parameters.data.numofgene - min(Truncations);
WrongPoints = nan(M,N);
SepCrit = nan(M,N);


%Construct Principal Eigenspace of Class A
NA = size(Datas.A.Training,2);
ZA = 1/sqrt(NA - 1)*Datas.A.Training;
[PhiA,SA] = mysvd(ZA, max(Truncations));

NB = size(Datas.B.Training,2);
ZB = 1/sqrt(NB - 1) * (Datas.B.Training - mean(Datas.B.Training,2));

parfor ima = Truncations
%for ima = Truncations

    UA = PhiA(:,1:ima);
    G = ZB * (ZB' - (ZB'*UA)*UA');
    Mres = P - ima;

    %Get a basis of eigenvectors corresponding to G

%     tic 
%     [T,~] = eig(G);
%     t1 = toc;
%     disp(toc)

    tic
    rankG = min([NB-1, Mres]);
    [T1,~] = eigs(G,rankG);
    G2 = T1' - (T1'*UA)*UA';
    T2 = null(G2);
    %norm(G*T2);
    %T3 = orth()
    T = [T1 , T2];
    t2 = toc;
    disp(t2);

    if ~isreal(T), keyboard, end

%     tic
%     [PhiA,~,~] = svd(ZA*ZA');
%     t3 = toc;
%     disp(t3)

    switch parameters.multilevel.eigentag
        case 'largest'
        case 'smallest'
            T = fliplr(T);
    end
    
    W = nan([1,N]);

    for imres = 10:Mres
        Si = T(:,1:imres);
        W(imres) = IdentifyMisplaced(Si, Datas, parameters);

    end

    WrongPoints(ima,:) = W;
end



end

function Datas = InScriptPrepData(Datas, parameters, methods)

%nTest = floor(min(parameters.data.A, parameters.data.B)/10);
nTest = parameters.Kfold;
if isempty(nTest)
    m = min(parameters.data.A, parameters.data.B);
    nTest = floor(m/10);
end

Datas.A.Training = Datas.rawdata.AData(:,1:end-nTest);
Datas.B.Training = Datas.rawdata.BData(:,1:end-nTest);
Datas.A.Testing = Datas.rawdata.AData(:,end-nTest+1:end);
Datas.B.Testing = Datas.rawdata.BData(:,end-nTest+1:end);

NA = size(Datas.A.Training,2);
NB = size(Datas.B.Training,2);


iData = 1:NA;
if parameters.multilevel.splitTraining
    %nTesting = NB;
    iCov = iData(iData <= NB);
else
    iCov = iData;
end

meanXA = mean(Datas.A.Training(:,iCov), 2);
for i = 'AB', for j = ["Training", "Testing"]
        Datas.(i).(j) = Datas.(i).(j) - meanXA;
end, end

Datas.A.Training = Datas.A.Training - meanXA;
Datas.B.Training = Datas.B.Training - meanXA ; 
Datas.A.Testing = Datas.A.Testing - meanXA;
Datas.B.Testing = Datas.B.Testing - meanXA;
end
% ========================================================================

%==========================================================================
function Truncations = GetTruncations(parameters)

if ~isempty(parameters.snapshots.k1)
    Truncations = parameters.snapshots.k1;
    return
end

%Get maximum truncation 
switch parameters.multilevel.splitTraining
    case true
        minTrainingA = parameters.data.A - parameters.Kfold;
        minTestingB = max(parameters.Kfold, mod(parameters.data.B, parameters.Kfold) );
        maxTrainingB = parameters.data.B - minTestingB;
        maxTrunc = minTrainingA - maxTrainingB;
        maxTrunc = min(parameters.data.numofgene-1, maxTrunc);
    case false
        maxTrunc = parameters.data.numofgene - 1;
end

%Get list of truncation parameters
Truncations = 1:maxTrunc;

end 
%==========================================================================

function wrong = IdentifyMisplaced(Si, Datas, parameters)

        XA = Si' * Datas.A.Training;
        XB = Si' * Datas.B.Training;
        mB = mean(XB,2);

        [uA, sA] = mysvd(XA);
        [uB, sB] = mysvd(XB);

        %KAinv = uA * ((sA.^(-0.5)) .* uA'); 
        %KBinv = uB * ((sB.^(-0.5)) .* uB'); 

        KAinv = ((sA.^(-0.5)) .* uA');
        KBinv = ((sB.^(-0.5)) .* uB'); 

        xA = KAinv * XA;
        xB = KBinv * (XB - mB);

        rA = sum(xA.^2,1);
        rB = sum(xB.^2,1);

        if ~isreal(rA), keyboard, end
        if ~isreal(rB), keyboard, end

        radA = quantile(rA, parameters.multilevel.concentration);
        radB = quantile(rB, parameters.multilevel.concentration);

        BinA = KAinv * XB;
        AinB = KAinv * (XA - mB);

        wrongA = sum(BinA.^2,1) < radA;
        wrongB = sum(AinB.^2,1) < radB;

        wrong = sum(wrongA) + sum(wrongB);

end