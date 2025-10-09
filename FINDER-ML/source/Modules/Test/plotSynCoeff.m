% plots coefficients from the synthetic data -- we expect the
% coefficients of the B class to be zero since the multilevel filter
% is generated from this class.
function plotSynCoeff(parameters, Datas)

levelcoeff = parameters.Coeff(1).B.levelcoeff;
maxlevel = max(levelcoeff);
numlevel = maxlevel;
counter = 1;


Q = Datas.B.Training(:,1);
k = parameters.data.nk;

for i=1:k
    figure(i)
    coeff = parameters.Coeff(i).B.C(:,1);

    subplot(maxlevel + 2,1,counter);
    stem(Q);

    for n = maxlevel : -1 : maxlevel - maxlevel
        counter = counter + 1;
        subplot(numlevel + 2, 1, counter);  
        stem(coeff(levelcoeff == n));
        title(['Level Coefficients = ',num2str(n)]);
    end
    
    sgtitle([num2str(parameters.snapshots.Brs(i)), 'points']);   
    counter=1;
end
end

