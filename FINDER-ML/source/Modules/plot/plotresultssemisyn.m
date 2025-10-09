function plotresultssemisyn(results, parameters)

   k = parameters.data.nk;
   j=parameters.multilevel.l;
   x=0:1:j;
   
   
   ml_acc=results.multilevel.accuracy;
   ml_prec=results.multilevel.precision;
   
   svm_acc=results.svm.accuracy;
   svm_prec=results.svm.precision;
   
   figure 
   tiledlayout(2,1)
   nexttile
   for i = 1:k
       plot(x,ml_acc(i,:), 'o-')
       hold on
   end
   hold off 
   legend('150', '250','450','850','1500');
   title('Multilevel Accuracy');
   xlabel('Level');
   nexttile
      for i = 1:k
       plot(x,ml_prec(i,:), 'o-')
       hold on
   end
   hold off 
   legend('150', '250','450','850','1500');
   title('Multilevel Precision');
   xlabel('Level');
   
   
   x=parameters.snapshots.normrs;
   figure 
   tiledlayout(2,1)
   nexttile
  
   plot(x,svm_acc, 'o-')

   title('SVM Accuracy');
   xlabel('# Simulated Points');
   nexttile
    
       plot(x,svm_prec, 'o-')
 
   title('SVM Precision');
   xlabel('# Simulated Points');
   

end 