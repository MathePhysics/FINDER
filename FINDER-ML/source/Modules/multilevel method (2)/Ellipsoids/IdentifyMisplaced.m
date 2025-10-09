%==========================================================================
function wrong = IdentifyMisplaced(Datas, parameters, methods)

        nargoutchk(1,3);

        %% Get info on the principal axes of each Class
        for C = 'AB' 
            NC = 1/ sqrt( size(Datas.(C).Machine,2) - 1);
            MC = mean(Datas.(C).Machine, 2);
            XC = NC * (Datas.(C).Machine - MC);
            [E.(C).UC, E.(C).SC] = svd(XC, 'econ', 'vector'); 
            %eigendata.(['Eval' C]) = E.(C).SC;
            %eigendata.(['Evec' C]) = E.(C).UC;
            E.(C).center = MC;
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
                E.(C).Kinv = (SC.^(-0.5) .* UC');
            elseif isempty(SC)
                E.(C).Kinv = 1;
            end
        end


        %% Estimate the "radius" of each class ellipse
        for C = 'AB'
            %XC = [Datas.(C).Machine, Datas.(C).Testing];
            %XC = Datas.(C).Testing;
            XC = Datas.(C).Machine;
            ZC = E.(C).Kinv*( XC - E.(C).center);
            radii = sum(ZC.^2, 1);
            E.(C).radius = quantile(radii, parameters.multilevel.concentration);
        end

        
        %% Find the number of wrongly misplaced points
        wrong = [];
        for Set = ["Machine", "Testing"]
            switch Set
                case "Machine", weight = 1; 
                case "Testing", weight = 3;
            end

        for C = 'AB'
            % XC = [Datas.(C).Machine, Datas.(C).Testing];
            for D = 'AB'
            %XC = Datas.(C).Testing;
            XC = Datas.(C).(Set);
            %switch C, case 'A', D = 'B'; case 'B', D = 'A'; end           
            ZC = E.(D).Kinv*(XC - E.(D).center);
            radii = sum(ZC.^2,1);
            switch C == D
                case false, iswrong = weight*4*sum(radii <= E.(D).radius);
                case true, iswrong = weight*sum(radii > E.(D).radius); 
            end
            wrong = [wrong , iswrong];
            
            end

        end
        end
        
        wrong = floor(sum(wrong));

end
        %wrong = max(wrong);
        
        %wrong = floor(sqrt(wrong));
%==========================================================================