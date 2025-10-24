function Datas = ProjectOntoT(Datas, parameters, methods)

%% Construct a Basis of R^(P-MA) using the multilevel tree method and project the Semi-Filtered Class A and B machine data onto this basis
    
% NA = size(Datas.A.Machine, 2);
% NB = size(Datas.B.CovTraining, 2);

NB = size(Datas.B.CovTraining,2);
%XB = 1/sqrt(NB - 1)*(Datas.B.CovTraining - mean(Datas.B.CovTraining,2));
XB = Datas.B.CovTraining;
XB = XB - mean(XB,2);

ProjectOnto = @(X,Y,Z,P) Coeffhbtrans(X, Y, Z, ...
                        P.Training.origin.multilevel.multileveltree, ...
                        P.Training.origin.multilevel.ind, ...
                        P.Training.origin.multilevel.datacell, ...
                        P.Training.origin.multilevel.datalevel);

m = size(XB);
    if ~methods.Multi2.isTallMatrix(XB) %m(2) >= 0.5*m(1)
        [T,~,~] = svd(XB);
        T = fliplr(T);
        %Troubleshooting
         %T = T';


        for C = 'AB', for set = ["Machine", "Testing"]
                Datas.(C).set = T'*Datas.(C).(set);
        end, end

    

    else
        p3 = parameters;
        iend = min([m(1), m(2) - 1]);
        p3.snapshots.k1 = min( iend );
        p3.Training.origin = methods.Multi.snapshots(XB, p3, methods, p3.snapshots.k1);
        p3 = methods.Multi.dsgnmatrix(methods, p3);
        p3 = methods.Multi.multilevel(methods, p3);

        %Troubleshooting
        % I = eye(size(T));
        % NC = size(I,2);
        % ZC = zeros(p3.Training.dsgnmatrix.origin.numofpoints, NC+1);
        % [I, LC, ~,~] = ProjectOnto(NC, ZC, I, p3);
        % I(:,end) = [];
        % I(LC == -1, :) =  p3.Training.origin.snapshots.eigenfunction(iend:-1:1,:);
        % 
        % I0 = I(LC == -1,:);
        % T0 = T(LC == -1,:);
        % Ip = I(LC > -1,:);
        % Tp = T(LC > -1,:);
        % 
        % disp(norm(I0' - T0' * (T0*I0'),'fro'));
        % disp(norm(Ip' - Tp' * (Tp*Ip'),'fro'));
        % disp(norm(T0' - I0' * (I0*T0'),'fro'));
        % disp(norm(Tp' - Ip' * (Ip*Tp'),'fro'));
        % I0T0diff = sum(abs(I0 - flipud(T0)),2);
        % plot(I0T0diff)
        
        
        for C = 'AB', for set = ["Machine", "Testing"]
          raw = Datas.(C).(set);
          NC = size(Datas.(C).(set),2);
          ZC = zeros(p3.Training.dsgnmatrix.origin.numofpoints, NC + 1);
          [Datas.(C).(set), LC, ~,~] = ProjectOnto(NC, ZC, Datas.(C).(set), p3);
          Datas.(C).(set)(:,end) = [];
          Datas.(C).(set)(LC == -1, :) =  p3.Training.origin.snapshots.eigenfunction(iend:-1:1,:) * raw;
          
        end, end


    end

end
 