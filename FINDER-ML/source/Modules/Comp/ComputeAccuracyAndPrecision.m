function results = ComputeAccuracyAndPrecision(results, l)

narginchk(1,2);

precision = (results.correct_A + results.correct_B)/...
    (results.correct_A + results.wrong_B + results.correct_B + results.wrong_A);
accuracy = results.correct_A/(results.correct_A + results.wrong_A);

if nargin == 2
    results.multilevel.precision(l) = precision;
    results.multilevel.accuracy(l) = accuracy;

elseif nargin == 1
    results.svm.precision = precision;
    results.svm.accuracy = accuracy;
end

end