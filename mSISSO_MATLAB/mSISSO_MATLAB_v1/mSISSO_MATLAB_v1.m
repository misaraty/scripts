clear;
clc;
close all
load('sissoTestData.mat')
fprintf("Fitting 'small' data: %d data points, %d features.\n", length(y), size(x, 2))
nNonzeroCoefs = 3;
nFeaturesPerSisIter = 10;
fprintf("Searching for models up to %d dimemsions, considering %d new features per iteration.\n", nNonzeroCoefs, nFeaturesPerSisIter)
sisso = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter, true);
sisso = fitSisso(sisso, x, y);
printModels(sisso, featureList)
load('sissoTestDataBig.mat')
fprintf("Fitting 'big' data: %d data points, %d features.\n", length(y), size(x, 2))
nNonzeroCoefs = 3;
nFeaturesPerSisIter = 26;
fprintf("Searching for models up to %d dimemsions, considering %d new features per iteration.\n", nNonzeroCoefs, nFeaturesPerSisIter)
sissoBig = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter, true);
sissoBig = fitSisso(sissoBig, x, y);
printModels(sissoBig, featureList)

function obj = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter, allL0Combinations)
if nargin < 1 || isempty(nNonzeroCoefs), nNonzeroCoefs = 1; end
if nargin < 2 || isempty(nFeaturesPerSisIter), nFeaturesPerSisIter = 1; end
if nargin < 3 || isempty(allL0Combinations), allL0Combinations = true; end
obj.nNonzeroCoefs = nNonzeroCoefs;
obj.nFeaturesPerSisIter = nFeaturesPerSisIter;
obj.allL0Combinations = allL0Combinations;
obj.coefs = [];
obj.coefsStandardized = [];
obj.intercept = [];
obj.listOfCoefs = {};
obj.listOfNonzeroCoefs = {};
obj.listOfIntercepts = {};
obj.rmses = [];
obj.selectedIndicesSis = {};
obj.unselectedIndicesSis = [];
obj.selectedIndicesL0 = {};
obj.selectedIndicesCurrent = [];
obj.scales = 1;
obj.means = 0;
end
function obj = fitSisso(obj, x, y)
warning('off')
checkParams(obj, size(x, 2))
obj.unselectedIndicesSis = 1:size(x, 2);
obj = setMeansAndScales(obj, x);
x = standardizeData(obj, x);
yMean = mean(y);
yCentered = y - yMean;
for iter = 1:obj.nNonzeroCoefs
    if iter == 1
        residuals = yCentered;
    else
        residuals = yCentered - predictFromStandardizedData(obj, x);
    end
    [obj, nClosestIndices, bestProjectionScore] = sis(obj, x, residuals);
    obj.selectedIndicesSis = [obj.selectedIndicesSis, {nClosestIndices}];
    if iter == 1
        obj.coefsStandardized = bestProjectionScore / length(y);
        obj.selectedIndicesCurrent = nClosestIndices(1);
        rmse = sqrt(sum((yCentered - predictFromStandardizedData(obj, x)).^2)./length(y));
    else
        [obj.coefsStandardized, obj.selectedIndicesCurrent, rmse] = exhaustiveRegularization(obj, x, yCentered, obj.selectedIndicesSis);
    end
    [coefsThisIter, obj.intercept] = getNonstandardizedCoefs(obj, obj.coefsStandardized, yMean);
    obj.coefs = zeros(1, size(x, 2));
    obj.coefs(obj.selectedIndicesCurrent) = coefsThisIter;
    obj.listOfNonzeroCoefs = [obj.listOfNonzeroCoefs, {coefsThisIter}];
    obj.listOfCoefs = [obj.listOfCoefs, {obj.coefs}];
    obj.listOfIntercepts = [obj.listOfIntercepts, {obj.intercept}];
    obj.selectedIndicesL0 = [obj.selectedIndicesL0, {obj.selectedIndicesCurrent}];
    obj.rmses = [obj.rmses, rmse];
end
end
function yPred = predictSisso(obj, x, dim)
if nargin < 3 || isempty(dim), dim = obj.nNonzeroCoefs; end
xModel = x(:, obj.selectedIndicesL0{dim});
xModel = [ones(size(x, 1), 1), xModel];
coefsModel = [obj.listOfIntercepts{dim}; obj.listOfNonzeroCoefs{dim}];
yPred = xModel * coefsModel;
end
function printModels(obj, features)
fprintf("%14s %16s\n", 'RMSE', 'Model')
for modelDim = 1:obj.nNonzeroCoefs
    disp(getModelString(obj, features, modelDim))
end
disp(" ")
end
function checkParams(obj, nColumns)
if nColumns < (obj.nNonzeroCoefs * obj.nFeaturesPerSisIter)
    error("nNonzeroCoefs*nFeaturesPerSisIter is larger than the number of columns in your input matrix. Choose smaller nNonzeroCoefs or nFeaturesPerSisIter.")
end
end
function obj = setMeansAndScales(obj, x)
obj.means = mean(x);
obj.scales = std(x);
obj.scales(obj.scales == 0) = 1;
end
function x = standardizeData(obj, x)
x = (x - obj.means) ./ obj.scales;
end
function yPred = predictFromStandardizedData(obj, x)
yPred = x(:, obj.selectedIndicesCurrent) * obj.coefsStandardized;
end
function [obj, indicesNClosestOut, bestProjectionScore] = sis(obj, x, y)
projectionScores = y' * x(:, obj.unselectedIndicesSis);
absProjectionScores = abs(projectionScores);
[~, indicesSorted] = sort(absProjectionScores, 'descend');
indicesNClosest = indicesSorted(1:obj.nFeaturesPerSisIter);
bestProjectionScore = projectionScores(indicesNClosest(1));
indicesNClosestOut = obj.unselectedIndicesSis(indicesNClosest);
mask = true(size(obj.unselectedIndicesSis));
mask(indicesNClosest) = false;
obj.unselectedIndicesSis = obj.unselectedIndicesSis(mask);
end
function [coefsStandardized, selectedIndicesCurrent, rmse] = exhaustiveRegularization(obj, x, y, listOfSisIndices)
squareErrorMin = y' * y;
coefsMin = [];
indicesCombinationMin = [];
if obj.allL0Combinations
    allIdx = [listOfSisIndices{:}];
    combinations = nchoosek(allIdx, length(listOfSisIndices));
else
    combinations = cartesianProduct(listOfSisIndices);
end
for i = 1:size(combinations, 1)
    thisCombination = combinations(i, :);
    xCombination = x(:, thisCombination);
    coefsCombination = xCombination \ y;
    residuals = y - (xCombination * coefsCombination);
    squareError = sum(residuals.^2);
    if squareError < squareErrorMin
        squareErrorMin = squareError;
        coefsMin = coefsCombination;
        indicesCombinationMin = thisCombination;
    end
end
coefsStandardized = coefsMin;
selectedIndicesCurrent = indicesCombinationMin;
rmse = sqrt(squareErrorMin/size(x, 1));
end
function [coefsOrig, interceptOrig] = getNonstandardizedCoefs(obj, coefsStandardized, interceptStandardized)
coefsOrig = coefsStandardized ./ obj.scales(obj.selectedIndicesCurrent)';
interceptOrig = interceptStandardized - (obj.means(obj.selectedIndicesCurrent) ./ obj.scales(obj.selectedIndicesCurrent)) * coefsStandardized;
end
function modelstr = getModelString(obj, features, modelDim)
coefs_dim = [obj.listOfIntercepts{modelDim}, obj.listOfNonzeroCoefs{modelDim}'];
selectedFeatures = {};
for k = obj.selectedIndicesL0{modelDim}
    selectedFeatures = [selectedFeatures, features(k)];
end
modelstr = sprintf("%dD:\t%8f\t", modelDim, obj.rmses(modelDim));
for dim = 1:modelDim + 1
    c = coefs_dim(dim);
    c_abs = abs(c);
    if dim == 1
        modelstr = strcat(modelstr, sprintf("%0.3f ", c));
    else
        if c >= 0, sgn = '+';
        else, sgn = '-';
        end
        modelstr = strcat(modelstr, sprintf("%s %0.3f %s ", sgn, c_abs, selectedFeatures{dim-1}));
    end
end
end
function x = cartesianProduct(sets)
numSets = length(sets);
sizeThisSet = zeros(numSets, 1);
for i = 1:numSets
    thisSet = sort(sets{i});
    if ~isvector(thisSet), error('All inputs must be vectors.'); end
    if ~isnumeric(thisSet), error('All inputs must be numeric.'); end
    if ~isequal(thisSet, unique(thisSet)), error(['Input set ', num2str(i), ' contains duplicated elements.']); end
    sizeThisSet(i) = length(thisSet);
    sets{i} = thisSet;
end
x = zeros(prod(sizeThisSet), numSets);
for i = 1:size(x, 1)
    idxVect = ind2subVect(sizeThisSet, i);
    for j = 1:numSets
        x(i, j) = sets{j}(idxVect(j));
    end
end
end
function X = ind2subVect(siz, idx)
if iscolumn(siz), siz = siz'; end
n = length(siz);
k = [1, cumprod(siz(1:end-1))];
idx = idx - 1;
X = zeros(1, n);
for i = n:-1:1
    X(i) = floor(idx/k(i)) + 1;
    idx = rem(idx, k(i));
end
end
function XOut = generateDescriptors(x, xVars)
for idxGroup = 1:length(x)
    [A1_X, A1_Xvars] = operatorA1(x{idxGroup}, xVars{idxGroup});
    [A2_X, A2_Xvars] = operatorA2(x{idxGroup}, xVars{idxGroup});
    [A3_X, A3_Xvars] = operatorA3(x{idxGroup}, xVars{idxGroup});
    [A4_X, A4_Xvars] = operatorA4(x{idxGroup}, xVars{idxGroup});
    [A5_X, A5_Xvars] = operatorA5(x{idxGroup}, xVars{idxGroup});
    [A6_X, A6_Xvars] = operatorA6(x{idxGroup}, xVars{idxGroup});
    x{idxGroup} = [x{idxGroup}, A1_X, A2_X, A3_X, A4_X, A5_X, A6_X];
    xVars{idxGroup} = [xVars{idxGroup}, A1_Xvars, A2_Xvars, A3_Xvars, A4_Xvars, A5_Xvars, A6_Xvars];
    [B1_X, B1_Xvars] = operatorB1(x{idxGroup}, xVars{idxGroup});
    x{idxGroup} = [x{idxGroup}, B1_X];
    xVars{idxGroup} = [xVars{idxGroup}, B1_Xvars];
end
for idxGroup = 1:length(x)
    if size(x{idxGroup}, 2) > 1
        [C2_X, C2_Xvars] = operatorC2(x{idxGroup}, xVars{idxGroup});
        x{idxGroup} = [x{idxGroup}, C2_X];
        xVars{idxGroup} = [xVars{idxGroup}, C2_Xvars];
    end
end
combinations = nchoosek(1:length(x), 2);
for idxCombination = 1:size(combinations, 1)
    thisCombination = combinations(idxCombination, :);
    [x{end+1}, xVars{end+1}] = operatorC1(x{thisCombination(1)}, x{thisCombination(2)}, xVars{thisCombination(1)}, xVars{thisCombination(2)});
end
x = [x{:}];
xVars = [xVars{:}];
mask = all(isfinite(x), 1);
Xout = x(:, mask);
Xoutvars = xVars(mask);
XOut = array2table(Xout, 'VariableNames', Xoutvars);
end
function [Xout, Xoutvars] = operatorA1(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^0.5;
    Xoutvars{i} = strcat('(', Xvars{i}, '^0.5', ')');
end
end
function [Xout, Xoutvars] = operatorA2(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^2;
    Xoutvars{i} = strcat('(', Xvars{i}, '^2', ')');
end
end
function [Xout, Xoutvars] = operatorA3(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^3;
    Xoutvars{i} = strcat('(', Xvars{i}, '^3', ')');
end
end
function [Xout, Xoutvars] = operatorA4(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^(1 / 3);
    Xoutvars{i} = strcat('(', Xvars{i}, '^(1/3)', ')');
end
end
function [Xout, Xoutvars] = operatorA5(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^(1 / 4);
    Xoutvars{i} = strcat('(', Xvars{i}, '^(1/4)', ')');
end
end
function [Xout, Xoutvars] = operatorA6(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = X(:, i).^4;
    Xoutvars{i} = strcat('(', Xvars{i}, '^4', ')');
end
end
function [Xout, Xoutvars] = operatorB1(X, Xvars)
Xout = X;
Xoutvars = Xvars;
for i = 1:size(X, 2)
    Xout(:, i) = 1 ./ X(:, i);
    Xoutvars{i} = strcat('(1/', Xvars{i}, ')');
end
end
function [Xout, Xoutvars] = operatorC1(X1, X2, X1vars, X2vars)
numdatapoints = size(X1, 1);
numfeaturevars = size(X1, 2) * size(X2, 2);
Xout = zeros(numdatapoints, numfeaturevars);
Xoutvars = cell(1, numfeaturevars);
idxXout = 1;
for i = 1:size(X1, 2)
    for j = 1:size(X2, 2)
        Xout(:, idxXout) = X1(:, i) .* X2(:, j);
        Xoutvars{idxXout} = strcat(X1vars{i}, '*', X2vars{j});
        idxXout = idxXout + 1;
    end
end
end
function [Xout, Xoutvars] = operatorC2(X1, X1vars)
numdatapoints = size(X1, 1);
pairs = nchoosek(1:size(X1, 2), 2);
numfeaturevars = size(pairs, 1);
Xout = zeros(numdatapoints, numfeaturevars);
Xoutvars = cell(1, numfeaturevars);
for i = 1:numfeaturevars
    a = pairs(i, 1);
    b = pairs(i, 2);
    Xout(:, i) = abs(X1(:, a)-X1(:, b));
    Xoutvars{i} = strcat('|', X1vars{a}, '-', X1vars{b}, '|');
end
end
