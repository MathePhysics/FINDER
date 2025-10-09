function plotcoeff(parameters, results, Datas);

levelcoeff = results.Training.normal.levelcoeff;
maxlevel = max(levelcoeff);
numlevel = maxlevel;
counter = 1;

x = 0 : ( 1 / (parameters.Training.dsgnmatrix.origin.numofpoints-1) ) : 1;
Q = Datas.normal.Training(:,1);

coeff = results.Training.normal.C(:,1);


subplot(maxlevel + 2,1,counter);
plot(x,Q);


for n = maxlevel : -1 : maxlevel - maxlevel
    counter = counter + 1;
    subplot(numlevel + 2, 1, counter);  
    stem(coeff(levelcoeff == n));
    title(['Level Coefficients = ',num2str(n)]);
end


end

