function [C, levelcoeff, dcoeffs, ccoeffs] = Coeffhbtrans(m, C, Q, multileveltree, inds, datacell, datalevel);

%parfor i = 1:m
for i = 1:m
    hQ = Q(:,i);
    [coeff, levelcoeff, dcoeffs, ccoeffs] = hbtrans(hQ, multileveltree, inds, datacell, datalevel);
    C(:,i) = coeff;
    
     
    
end

C(:,m+1) = levelcoeff;


end