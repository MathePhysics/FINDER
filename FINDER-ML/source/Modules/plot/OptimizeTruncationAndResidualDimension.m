function OptimizeTruncationAndResidualDimension(Datas, parameters, methods, results)

map.Sep = zeros(parameters.data.numofgene); 
map.Close = map.Sep;

Datas.A.Training = Datas.rawdata.AData;
Datas.B.Training = Datas.rawdata.BData;
Datas.A.Testing = sparse(parameters.data.numofgene, 1);
Datas.B.Testing = Datas.A.Testing;

ZA = Datas.A.Training - mean(Datas.A.Training,2);
[Datas.A.eigenvectors, ~, ~] = svd(ZA);

for ima = 1:parameters.data.numofgene
    for imres = 1:ima
        results.current_ima = ima; results.current_imres = imres;



        for machine = ["Sep", "Close"]

        Datas.A.Training = Datas.rawdata.AData;
        Datas.B.Training = Datas.rawdata.BData;
        Datas.A.Testing = sparse(parameters.data.numofgene, 1);
        Datas.B.Testing = Datas.A.Testing;

        Datas = ConstructOptimalBasis(Datas,parameters,methods,results,machine);


        r = norm(mean(Datas.B.Training, 2)); %= distance between class means
        numA = mean(sum(Datas.A.Training.^2,1));
        ZB = Datas.B.Training - mean(Datas.B.Training, 2);
        numB = mean( sum(ZB.^2, 1) );
        ErrorBoundFunction = @(rA) numA / rA^2 + numB / (r - rA)^2 ;
        switch machine
            case "Sep"
            ErrorBound = fminbnd(ErrorBoundFunction,0,r);
            
            
            case "Close"
            ErrorBound = numB / numA;
            %fprintf('Ratio of Class B concentration to Class A: %0.3e \n\n', ErrorBound);
        end

        map.(machine)(ima, imres) = ErrorBound;



        end
    end
end

%[bestEB, ibest] = min(map.Sep( map.Sep ~= 0));
[bestEB, ibest] = max(map.Sep( map.Sep ~= 0));
[MA, Mres] = ind2sub(size(map.Sep), ibest);
fprintf('Separate means assumed \n Type I + Type II error = %0.3e \n', bestEB);
fprintf('Best MA: %d, Best Mres: %d \n\n', MA, Mres)

[bestRatio, ibest] = max(map.Close(:));
[MA, Mres] = ind2sub(size(map.Close), ibest);
fprintf('Close means assumed \n Type I : Type II ratio = %0.3e \n', bestRatio);
fprintf('Best MA: %d, Best Mres: %d \n\n', MA, Mres)

figure()
machines = {'Sep', 'Close'};

J = jet; J(1,:) = [0 0 0];
minmax = @(x) [min(x), max(x)];
for i = 1:2
    machine = machines{i};
    subplot(1,2,i)
    imagesc(map.(machine))
    c(:,i) = minmax(map.(machine));
    xlabel('M_{res}', 'Interpreter', 'tex'), ylabel('M_A', 'Interpreter', 'tex')
    title(machine)
    colorbar 
    colormap(J)
end

% c2 = [min(c(:)), max(c(:))];
% for i = 1:2
% set(subplot(1,2,i), 'clim', c2)
% end

