%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KLprojection
% Obtain residuals after projection onto space spanned by 
% first KL components and mean
% 
%% Arguments:
% -`parameters` and `data` structures from
% -`methods.domaindata` and `methods.KLeigenspace`
%
%% Value:  
% "parameters.KL.Output.residual" structure to parameters 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function parameters = KLprojection(methods,parameters,data)

%% input parameters ---------------------------------------
testdata_frame_number = parameters.data.testdata;
modality = parameters.data.modality;
M = parameters.KL.M; % eigenspace with mean

include_mean = parameters.operations.include_mean;
E = parameters.KL.mean;

%% transform testdata to vector ---------------------------
if isfield(data.landsat,'rawdata') 
    % older version
    testdata = reshape(data.landsat.rawdata(:,:,testdata_frame_number,modality), [], 1);
else
    % newer version (under testing)
    testdata = squeeze(data.landsat.reshapedata(:,testdata_frame_number,1));
end
% linear_indices_of_testdata = sub2ind(size(testdata),data.landsat.coords(:,1),data.landsat.coords(:,2));
linear_indices_of_testdata = data.landsat.pos; 
testdata_vector = testdata(linear_indices_of_testdata);


%% compute residual ---------------------------------------
if ~include_mean
    testdata_vector = testdata_vector - E;
end
residual = testdata_vector - M*(M'*testdata_vector); %parentheses speed computation time
testdata(linear_indices_of_testdata) = residual; 


%% save results in parameters structure -------------------
parameters.KL.Output.residual = testdata;

end