function totalerror = Testchange(coeff, dcoeffs, ccoeffs, multileveltree, ind, datacell, datalevel, numofpoints);


% Test format change from vector of coefficients to struct
totalerror = 0;

[a, b] = hbvectortocoeffs(coeff, multileveltree, ind, datacell, datalevel, numofpoints);

for i = 1 : length(dcoeffs)
   totalerror = totalerror + (norm(a{i} -  dcoeffs{i}));
end

totalerror = totalerror + (norm(b -  ccoeffs));
fprintf('Total error = %e \n', totalerror);

end




