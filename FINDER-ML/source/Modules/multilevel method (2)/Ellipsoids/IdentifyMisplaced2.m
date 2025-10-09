%==========================================================================
function wrong = IdentifyMisplaced2(Datas, parameters, methods)

        nargoutchk(1,3);

        %% Get info on the principal axes of each Class
        for C = 'AB' 
            NC = 1/ sqrt( size(Datas.(C).Machine,2) - 1);
            MC = mean(Datas.(C).Machine, 2);
            XC = NC * (Datas.(C).Machine - MC);
            [E.(C).UC, E.(C).SC, ~] = svd(XC, 'econ', 'vector'); 
            E.(C).totvar = sum(E.(C).SC);
            %eigendata.(['Eval' C]) = E.(C).SC;
            %eigendata.(['Evec' C]) = E.(C).UC;
            E.(C).center = MC;
        end 

        %% Identify class with larger variance
        totvar = [E.A.totvar, E.B.totvar];
        [~,ic] = max(totvar);
        switch ic
            case 1, C1 = 'A'; C2 = 'B'; 
            case 2, C1 = 'B'; C2 = 'A';
        end

        %% Drop axes corresponding to sufficiently small singular values
        LSV = max([E.A.SC; E.B.SC]); %Largest Singular Value
        scaleFactor = 1/sqrt(max([size(E.A.UC,2), size(E.B.UC,2)]));
        zeroThresh = LSV * scaleFactor * eps;
        for C = 'AB'
            isSuffLarge = E.(C).SC >= zeroThresh;
            SC = E.(C).SC(isSuffLarge);
            UC = E.(C).UC(:,isSuffLarge);
            % Construct pseudo-inverse, which transforms data into
            % isotropic data
            if ~isempty(SC)
                %E.(C).Kinv = (SC.^(-0.5) .* UC');
                E.(C).Kinv = (SC.^-1) .* UC';
            elseif isempty(SC)
                E.(C).Kinv = 1;
            end
        end


        %% Estimate the "radius" of each class ellipse
        for C = 'AB'
            XC = [Datas.(C).Machine, Datas.(C).Testing];
            ZC = E.(C).Kinv*( XC - E.(C).center);
            radii = sum(ZC.^2, 1);
            E.(C).radius = quantile(radii, parameters.multilevel.concentration);
        end

        
        %% Find the number of wrongly misplaced points
        %A point will be considered misclassified if it is a member of the
        %class with larger spread, and lies in the other class's ellipse

        X = [Datas.(C1).Machine, Datas.(C1).Testing];
        Z = E.(C2).Kinv*(X - E.(C2).center);
        radii = sum(Z.^2,1);
        wrong = sum(radii < E.(C2).radius);
       



end
%==========================================================================