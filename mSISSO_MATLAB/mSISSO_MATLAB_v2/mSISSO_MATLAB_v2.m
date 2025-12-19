clear;
clc;
close all;
opts = detectImportOptions('data_v4.csv', 'VariableNamingRule', 'preserve');
T = readtable('data_v4.csv', opts);
X0 = T{:, 6:13};
% X0 = T{:, 6:9};
% X0 = T{:, 10:13};
% y = T{:, end};
y = T{:, end-1};
featureList0 = string(T.Properties.VariableNames(6:13));
fprintf("Loaded data_v4.csv: %d samples, %d base features.\n", size(X0, 1), size(X0, 2));
[X, featureList] = generateDescriptorLibrary(X0, featureList0);
fprintf("Generated descriptor library: %d samples, %d features.\n", size(X, 1), size(X, 2));
rng(42);
N = size(X, 1);
idx = randperm(N);
nTrain = round(0.70*N);
nVal = round(0.15*N);
nTest = N - nTrain - nVal;
idxTrain = idx(1:nTrain);
idxVal = idx(nTrain+1:nTrain+nVal);
idxTest = idx(nTrain+nVal+1:end);
Xtr = X(idxTrain, :);
ytr = y(idxTrain);
Xva = X(idxVal, :);
yva = y(idxVal);
Xte = X(idxTest, :);
yte = y(idxTest);
fprintf("Split: train=%d, val=%d, test=%d\n", nTrain, nVal, nTest);
nNonzeroCoefs = 3;
nFeaturesPerSisIter = 26;
fprintf("Searching for models up to %d dimensions, considering %d new features per iteration.\n", nNonzeroCoefs, nFeaturesPerSisIter);
sisso = SissoRegressor(nNonzeroCoefs, nFeaturesPerSisIter, true);
sisso = fitSisso(sisso, Xtr, ytr);
printModels(sisso, featureList);
dim = nNonzeroCoefs;
yhat_tr = predictSisso(sisso, Xtr, dim);
yhat_va = predictSisso(sisso, Xva, dim);
yhat_te = predictSisso(sisso, Xte, dim);
[mae_tr, rmse_tr, r2_tr] = regressionMetrics(ytr, yhat_tr);
[mae_va, rmse_va, r2_va] = regressionMetrics(yva, yhat_va);
[mae_te, rmse_te, r2_te] = regressionMetrics(yte, yhat_te);
fprintf("\n=== Metrics (Model Dim = %d) ===\n", dim);
fprintf("Train: MAE=%8.4f, RMSE=%8.4f, R2=%8.4f\n", mae_tr, rmse_tr, r2_tr);
fprintf("Val  : MAE=%8.4f, RMSE=%8.4f, R2=%8.4f\n", mae_va, rmse_va, r2_va);
fprintf("Test : MAE=%8.4f, RMSE=%8.4f, R2=%8.4f\n", mae_te, rmse_te, r2_te);
function [mae, rmse, r2] = regressionMetrics(yTrue, yPred)
yTrue = yTrue(:);
yPred = yPred(:);
err = yPred - yTrue;
mae = mean(abs(err));
rmse = sqrt(mean(err.^2));
ss_res = sum((yTrue - yPred).^2);
ss_tot = sum((yTrue - mean(yTrue)).^2);
if ss_tot == 0
    r2 = NaN;
else
    r2 = 1 - ss_res / ss_tot;
end
end
function [Xlib, namesLib] = generateDescriptorLibrary(X0, names0)
warning('off');
X0 = double(X0);
names0 = string(names0);
N = size(X0, 1);
eps0 = 1e-12;
Xun = [];
Nun = [];
Xun = [Xun, X0];
Nun = [Nun, names0];
Xun = [Xun, -X0];
Nun = [Nun, "(-" + names0 + ")"];
Xun = [Xun, abs(X0)];
Nun = [Nun, "abs(" + names0 + ")"];
Xun = [Xun, safePow(X0, 0.5)];
Nun = [Nun, "(" + names0 + "^0.5)"];
Xun = [Xun, safePow(X0, 2)];
Nun = [Nun, "(" + names0 + "^2)"];
Xun = [Xun, safePow(X0, 3)];
Nun = [Nun, "(" + names0 + "^3)"];
Xun = [Xun, safePow(X0, 4)];
Nun = [Nun, "(" + names0 + "^4)"];
Xun = [Xun, safePow(X0, 1/3)];
Nun = [Nun, "(" + names0 + "^(1/3))"];
Xun = [Xun, safePow(X0, 1/4)];
Nun = [Nun, "(" + names0 + "^(1/4))"];
Xun = [Xun, 1 ./ (X0 + eps0)];
Nun = [Nun, "(1/(" + names0 + "+eps))"];
Xun = [Xun, 1 ./ (abs(X0) + eps0)];
Nun = [Nun, "(1/(abs(" + names0 + ")+eps))"];
Xun = [Xun, safeLog(X0, eps0)];
Nun = [Nun, "log(" + names0 + "+eps)"];
Xun = [Xun, safeLog(abs(X0), eps0)];
Nun = [Nun, "log(abs(" + names0 + ")+eps)"];
Xun = [Xun, safeExp(X0)];
Nun = [Nun, "exp(" + names0 + ")"];
tmp = safeLog(abs(X0), eps0);
Xun = [Xun, tmp.^2];
Nun = [Nun, "(log(abs(" + names0 + ")+eps))^2"];
M = size(Xun, 2);
pairs = nchoosek(1:M, 2);
K = size(pairs, 1);
Xp = zeros(N, K);
Np = strings(1, K);
Xs = zeros(N, K);
Ns = strings(1, K);
Xsub = zeros(N, K);
Nsub = strings(1, K);
Xd = zeros(N, K);
Nd = strings(1, K);
Xr1 = zeros(N, K);
Nr1 = strings(1, K);
Xr2 = zeros(N, K);
Nr2 = strings(1, K);
Xnd = zeros(N, K);
Nnd = strings(1, K);
for k = 1:K
    i = pairs(k, 1);
    j = pairs(k, 2);
    a = Xun(:, i);
    b = Xun(:, j);
    Xp(:, k) = a .* b;
    Np(k) = Nun(i) + "*" + Nun(j);
    Xs(:, k) = a + b;
    Ns(k) = "(" + Nun(i) + "+" + Nun(j) + ")";
    Xsub(:, k) = a - b;
    Nsub(k) = "(" + Nun(i) + "-" + Nun(j) + ")";
    Xd(:, k) = abs(a-b);
    Nd(k) = "|" + Nun(i) + "-" + Nun(j) + "|";
    Xr1(:, k) = safeDiv(a, b, eps0);
    Nr1(k) = "(" + Nun(i) + "/(" + Nun(j) + "+eps))";
    Xr2(:, k) = safeDiv(b, a, eps0);
    Nr2(k) = "(" + Nun(j) + "/(" + Nun(i) + "+eps))";
    Xnd(:, k) = (a - b) ./ (abs(a) + abs(b) + eps0);
    Nnd(k) = "((" + Nun(i) + "-" + Nun(j) + ")/(abs(" + Nun(i) + ")+abs(" + Nun(j) + ")+eps))";
end
M3 = min(size(Xun, 2), 24);
trip = nchoosek(1:M3, 3);
Kt = size(trip, 1);
Xt_sum = zeros(N, Kt);
Nt_sum = strings(1, Kt);
Xt_prod = zeros(N, Kt);
Nt_prod = strings(1, Kt);
Xt_ratio = zeros(N, Kt);
Nt_ratio = strings(1, Kt);
for k = 1:Kt
    i = trip(k, 1);
    j = trip(k, 2);
    m = trip(k, 3);
    a = Xun(:, i);
    b = Xun(:, j);
    c = Xun(:, m);
    Xt_sum(:, k) = a + b + c;
    Nt_sum(k) = "(" + Nun(i) + "+" + Nun(j) + "+" + Nun(m) + ")";
    Xt_prod(:, k) = a .* b .* c;
    Nt_prod(k) = "(" + Nun(i) + "*" + Nun(j) + "*" + Nun(m) + ")";
    Xt_ratio(:, k) = safeDiv(a.*b, c, eps0);
    Nt_ratio(k) = "((" + Nun(i) + "*" + Nun(j) + ")/(" + Nun(m) + "+eps))";
end
Xall = [Xun, Xp, Xs, Xsub, Xd, Xr1, Xr2, Xnd, Xt_sum, Xt_prod, Xt_ratio];
Nall = [Nun, Np, Ns, Nsub, Nd, Nr1, Nr2, Nnd, Nt_sum, Nt_prod, Nt_ratio];
finiteMask = all(isfinite(Xall), 1);
Xall = Xall(:, finiteMask);
Nall = Nall(finiteMask);
s = std(Xall, 0, 1);
varMask = s > 0;
Xall = Xall(:, varMask);
Nall = Nall(varMask);
[Xlib, namesLib] = dedupColumnsFast(Xall, Nall);
warning('on');
end
function Xp = safePow(X, p)
Xp = NaN(size(X));
if abs(p-round(p)) < 1e-12
    Xp = X.^p;
else
    mask = X >= 0;
    Xp(mask) = X(mask).^p;
end
end
function Xd = safeDiv(a, b, eps0)
Xd = a ./ (b + eps0);
end
function Xl = safeLog(X, eps0)
Xl = NaN(size(X));
mask = X >= 0;
Xl(mask) = log(X(mask)+eps0);
end
function Xe = safeExp(X)
Xc = max(min(X, 50), -50);
Xe = exp(Xc);
end
function [Xout, Nout] = dedupColumnsFast(X, N)
mu = mean(X, 1);
sd = std(X, 0, 1);
mn = min(X, [], 1);
mx = max(X, [], 1);
fp = round([mu; sd; mn; mx]'*1e10) / 1e10;
[~, ia] = unique(fp, 'rows', 'stable');
Xout = X(:, ia);
Nout = N(ia);
end
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
warning('on')
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
