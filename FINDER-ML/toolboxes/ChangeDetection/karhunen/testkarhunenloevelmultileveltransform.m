%% Data and parameters
clc;
clear all;
close all;

degree = []; % Dummy variable, not used in this code.
numofpoints = 100;
data = linspace(0,1,numofpoints)';
n = 10;
x = 0 : ( 1 / (numofpoints-1) ) : 1;
M = ones(size(x'));
Lp = 1;

%eigenfunctions;
for i =  2 : n
    if floor(i/2) == i/2 
        phi = (1 / i) * sin ( floor(i/2) * pi * x / Lp );
        %disp('even');
    else
        %disp('odd');
        phi = (1 / i) * cos ( floor(i/2) * pi * x / Lp ); 
    end
    M = [M phi'];
end
% Optional design matrix
polymodel.M = M;

% Number of columns in the design matrix
params.indexsetsize = n;

%%

% Create Multilevel Binary tree
fprintf('\n');
fprintf('Create KDd tree ---------------------------------\n');
tic;
[datatree, sortdata] = make_tree(data,@split_KD,params);
toc;

% Create multilevel basis
fprintf('\n');
fprintf('Create multilevel basis ------------------------\n');
tic;
[multileveltree, ind, datacell, datalevel]  = multilevelbasis(datatree, sortdata, degree, polymodel);
toc;

%% Run transform with random data
%Realization;
fprintf('\n');
fprintf('Test transform ---------------------------------\n');
Q = polymodel.M * rand(n,1);

figure(1);
plot(x,Q);
title('KL Realization');


%% Metrics
numvecs = size(Q,2);
maxerror = inf;
maxwaverror = inf;

tic;
for n = 1 : numvecs
    [coeff, levelcoeff, dcoeffs, ccoeffs] = hbtrans(Q(:,n), multileveltree, ind, datacell, datalevel);
    Qv = invhbtrans(dcoeffs, ccoeffs, multileveltree, ind, datacell, datalevel, numofpoints);
    polyerror(n) = norm(Q(:,n) - Qv, 2) / norm(Q(:,n));
    wavecoeffnorm(n) = norm(coeff(1 : end - params.indexsetsize),'inf')/norm(Q(:,n));
end
t = toc;

fprintf('\n');
fprintf('Results: Num of tests = %d, Max Relative Error = %e \n', numvecs, max(polyerror));
fprintf('         Relative HBcoefficientsnorm =  %e \n', max(wavecoeffnorm));
fprintf('Total timing = %d seconds \n',t);
tic;
[coeff, levelcoeff, dcoeffs, ccoeffs] = hbtrans(Q(:,n), multileveltree, ind, datacell, datalevel);
t = toc;
fprintf('One Hierarchical Basis transform time = %d seconds \n',t);

% Test format change from vector of coefficients to struct
totalerror = 0;
[a b] = hbvectortocoeffs(coeff, multileveltree, ind, datacell, datalevel, numofpoints);
for i = 1 : length(dcoeffs)
   totalerror = totalerror + (norm(a{i} -  dcoeffs{i}));
end
totalerror = totalerror + (norm(b -  ccoeffs));
fprintf('Total error = %e \n', totalerror);

%% Plot coefficients

figure(2)
maxlevel = max(levelcoeff);
numlevel = 3;
counter = 1;
subplot(numlevel + 2,1,counter);
plot(x,Q);
for n = maxlevel : -1 : maxlevel - numlevel
    counter = counter + 1;
    subplot(numlevel + 2, 1, counter);  
    stem(coeff(levelcoeff == n));
    title(['Level Coefficients = ',num2str(n)]);
end

%% Run trans
fprintf('\n');
fprintf('Add Bump to KL  --------------------------------\n');

% Random realization
% Q = polymodel.M * rand(size(polymodel.M,2),1);

% Add Gaussian "bump" to the data
Qkl = Q;
sigma = 0.0001;
maxbump = 0.01;
Q = Q + maxbump * exp( - ((x' - 0.65).^2)/sigma);

numvecs = size(Q,2);
maxerror = inf;
maxwaverror = inf;

tic;
for n = 1 : numvecs
    [coeff, levelcoeff, dcoeffs, ccoeffs] = hbtrans(Q(:,n), multileveltree, ind, datacell, datalevel);
    Qv = invhbtrans(dcoeffs, ccoeffs, multileveltree, ind, datacell, datalevel, numofpoints);
    polyerror(n) = norm(Q(:,n) - Qv, 2) / norm(Q(:,n));
    wavecoeffnorm(n) = norm(coeff(1 : end - params.indexsetsize),'inf')/norm(Q(:,n));
end
t = toc;

%%
figure(3);
subplot(3,1,1);
plot(x,Qkl);
axis([0 1 0 2])
title('KL Realization with Gaussian bump');
subplot(3,1,2);
plot(x, maxbump * exp( - ((x' - 0.65).^2)/sigma));
title('Gaussian bump');
subplot(3,1,3);
plot(x,Q);
axis([0 1 0 2])
title('KL Realization with Gaussian bump');

%%
figure(4)
maxlevel = max(levelcoeff);
numlevel = 3;
counter = 1;
subplot(numlevel + 2,1,counter);
plot(x, maxbump * exp( - ((x' - 0.65).^2)/sigma));
for n = maxlevel : -1 : maxlevel - numlevel
    counter = counter + 1;
    subplot(numlevel + 2, 1, counter);  
    stem(coeff(levelcoeff == n));
    title(['Level Coefficients = ',num2str(n)]);
end


