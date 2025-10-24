function [Datas] = TransformTree(Datas, parameters, methods)



if ~parameters.transform.ComputeTransform, return, end
%if parameters.transform.istransformed, return, end



%Construct Eigendata
eigendata = methods.transform.eigendata(Datas,parameters,'Training');
% eigendata_test = methods.transform.eigendata(Datas,parameters,'Testing');
% %imagesc(c), title('Testing Before'), colorbar, colormap jet, set(gca, 'ColorScale', 'log')
% 
% subplot(2,2,1)
% X = 1:eigendata.RankA;, Y=  1:eigendata.RankB;
% [c,~,~,~] = Optimize_Closeness_Ratio(eigendata);
% imagesc(c), title('Training Before'), colorbar, colormap jet, set(gca, 'ColorScale', 'log')
% 
% subplot(2,2,2)
% X = 1:eigendata_test.RankA;, Y =  1:eigendata_test.RankB;
% [c,~,~,~] = Optimize_Closeness_Ratio(eigendata_test);
% imagesc(c), title('Testing Before'), colorbar, colormap jet, set(gca, 'ColorScale', 'log')


%Reduce dimensionality
isHighRank = eigendata.RankA + eigendata.RankB > parameters.data.numofgene;



if isHighRank
    DecaysSlowly = methods.transform.decayRate(eigendata, parameters);
    if DecaysSlowly
    end

    eigendata = methods.transform.dimReduction(eigendata, parameters);
    
%Datas.A.Training = eigendata.EvecA * eigendata.EvecA' * Datas.A.Training;
%Datas.A.Testing = eigendata.EvecA * eigendata.EvecA' * Datas.A.Training;
end


%Construct Transformation
iOverlap = methods.transform.overlap(eigendata, parameters);

if ~any(iOverlap)
    Datas = methods.transform.LI(eigendata, Datas);
else
    Datas = methods.transform.LD(eigendata, Datas);
end

H = Datas.H;
E = Datas.E;
% 
Datas.B.Training = E * H' * Datas.B.Training;
Datas.B.Testing = E * H' * Datas.B.Testing;
Datas.A.Training = E * H' * Datas.A.Training;
Datas.A.Testing = E * H' * Datas.A.Testing;

%parameters.transform.istransformed = true;


%Construct Eigendata
% eigendata = methods.transform.eigendata(Datas,parameters,'Training');
% eigendata_test = methods.transform.eigendata(Datas,parameters,'Testing');
% 
% subplot(2,2,3)
% X = 1:eigendata.RankA;, Y=  1:eigendata.RankB;
% [c,~,~,~] = Optimize_Closeness_Ratio(eigendata);
% imagesc(X,Y,c), title('Training After'), colorbar, colormap jet, set(gca, 'Colorscale', 'log')
% 
% subplot(2,2,4)
% X = 1:eigendata_test.RankA;, Y =  1:eigendata_test.RankB;
% [c,~,~,~] = Optimize_Closeness_Ratio(eigendata_test);
% imagesc(c), title('Testing After'), colorbar, colormap jet, set(gca, 'Colorscale', 'log')
% 
% saveas(gcf, fullfile(parameters.datafolder, 'Separation Criterion.png'))

0;




%[Datas, parameters] = ApplyTransformation(Datas, parameters);


