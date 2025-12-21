rng(42)

scriptFull = mfilename("fullpath");
if strlength(scriptFull) == 0
    baseDir = pwd;
else
    baseDir = fileparts(scriptFull);
end

logFile = fullfile(baseDir, "run.log");
fid = fopen(logFile, "w");
if fid > 0
    fclose(fid);
end
diary(logFile); diary on

tStart = tic;

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

dataFile = fullfile(baseDir, "data_v4.csv");
if ~exist(dataFile, "file")
    dataFile = "data_v4.csv";
end

T = readtable(dataFile, VariableNamingRule="preserve");

X = T{:, 6:13};
Y = T{:, end};
X = double(X);
Y = double(Y);

featNames = string(T.Properties.VariableNames(6:13));

N = size(X, 1);
K = 5;
cv = cvpartition(N, "KFold", K);

valRatio = 0.2;

models = ["xgboost_style","lightgbm_style","random_forest"];

tuneCfg.nTrials = 60;
tuneCfg.patience = 18;
tuneCfg.minImprove = 1e-4;
tuneCfg.useTrainR2Constraint = true;
tuneCfg.maxTrainR2Drop = 0.10;
tuneCfg.trainSubsampleMax = 1600;

tuneCfg.gapPenalty = 0.35;
tuneCfg.trainPenalty = 0.08;
tuneCfg.r2Weight = 0.14;

fprintf("BaseDir: %s\n", baseDir)
fprintf("Data: %s | N=%d, D=%d | K=%d | valRatio=%.2f\n", dataFile, N, size(X,2), K, valRatio)
fprintf("nTrials=%d | patience=%d | gapPenalty=%.3f | r2Weight=%.3f\n", tuneCfg.nTrials, tuneCfg.patience, tuneCfg.gapPenalty, tuneCfg.r2Weight)

for mi = 1:numel(models)
    modelName = models(mi);

    fprintf("\n===========================\n");
    fprintf("Running model: %s\n", modelName);
    fprintf("===========================\n");

    R2_train = zeros(K, 1);
    R2_val = zeros(K, 1);
    R2_test = zeros(K, 1);
    MAE_train = zeros(K, 1);
    MAE_val = zeros(K, 1);
    MAE_test = zeros(K, 1);
    RMSE_train = zeros(K, 1);
    RMSE_val = zeros(K, 1);
    RMSE_test = zeros(K, 1);

    YTrain_all = cell(K, 1);
    YHatTrain_all = cell(K, 1);
    YVal_all = cell(K, 1);
    YHatVal_all = cell(K, 1);
    YTest_all = cell(K, 1);
    YHatTest_all = cell(K, 1);

    imp_all = zeros(K, size(X,2));

    foldModels = cell(K, 1);
    foldInfos = cell(K, 1);

    totalSteps = K;
    hwb = waitbar(0, "Starting...", "Name", "progress");
    cleanupObj = onCleanup(@()safeCloseWaitbar(hwb));

    sharedInfo = [];
    sharedInfoLogged = false;

    for k = 1:K
        fprintf("Fold %d / %d\n", k, K)

        idxTrainFull = find(training(cv, k));
        idxTest = find(test(cv, k));

        p = randperm(numel(idxTrainFull));
        idxTrainFull = idxTrainFull(p);

        nTrainFull = numel(idxTrainFull);
        nVal = max(1, round(valRatio*nTrainFull));

        idxVal = idxTrainFull(end-nVal+1:end);
        idxTrain = idxTrainFull(1:end-nVal);

        XTrain = X(idxTrain, :);
        YTrain = Y(idxTrain);
        XVal = X(idxVal, :);
        YVal = Y(idxVal);
        XTest = X(idxTest, :);
        YTest = Y(idxTest);

        muX = mean(XTrain, 1);
        sigX = std(XTrain, 0, 1);
        sigX(sigX == 0) = 1;

        XTrainN = (XTrain - muX) ./ sigX;
        XValN = (XVal - muX) ./ sigX;
        XTestN = (XTest - muX) ./ sigX;

        muY = mean(YTrain);
        sigY = std(YTrain);
        if sigY == 0
            sigY = 1;
        end

        YTrainN = (YTrain - muY) ./ sigY;
        YValN = (YVal - muY) ./ sigY;
        YTestN = (YTest - muY) ./ sigY;

        if isvalid(hwb)
            waitbar(min(1, (k - 1)/totalSteps), hwb, sprintf("Fold %d/%d | %s", k, K, modelName));
            drawnow limitrate
        end

        rng(2025 + 100000*mi + 1000*k)

        d = size(XTrainN, 2);

        nTr = size(XTrainN, 1);
        if nTr > tuneCfg.trainSubsampleMax
            idxSub = randperm(nTr, tuneCfg.trainSubsampleMax);
        else
            idxSub = 1:nTr;
        end
        XTrainSub = XTrainN(idxSub, :);
        YTrainSub = YTrainN(idxSub);

        if k == 1
            rng(777 + 10*mi)
            bestScore = inf;
            bestMdl = [];
            bestInfo = struct();
            bestTrainR2N = -inf;
            noImprove = 0;

            for tr = 1:tuneCfg.nTrials
                [mdl, info, yhatValN, yhatTrSubN] = train_one_trial(modelName, XTrainN, YTrainN, XValN, YValN, XTrainSub, YTrainSub, d, mi, k, tr);

                [trainR2N, trainRMSE] = r2_rmse(YTrainSub, yhatTrSubN);
                [valR2N, valRMSE] = r2_rmse(YValN, yhatValN);

                ok = true;
                if tuneCfg.useTrainR2Constraint
                    if bestTrainR2N > -inf
                        if trainR2N < bestTrainR2N - tuneCfg.maxTrainR2Drop
                            ok = false;
                        end
                    end
                end

                gap = max(0, trainRMSE - valRMSE) + max(0, valRMSE - trainRMSE);
                score = valRMSE + tuneCfg.gapPenalty*gap + tuneCfg.trainPenalty*max(0, (1-trainR2N)) - tuneCfg.r2Weight*valR2N;

                if ok && score < bestScore - tuneCfg.minImprove
                    bestScore = score;
                    bestMdl = mdl;
                    bestInfo = info;
                    if tuneCfg.useTrainR2Constraint
                        bestTrainR2N = max(bestTrainR2N, trainR2N);
                    end
                    noImprove = 0;
                else
                    noImprove = noImprove + 1;
                end

                if noImprove >= tuneCfg.patience
                    break
                end
            end

            if isempty(bestMdl)
                error("No valid model selected on fold %d for %s", k, modelName)
            end

            sharedInfo = bestInfo;
        end

        if isempty(sharedInfo)
            error("Shared hyperparameters missing for %s", modelName)
        end

        rng(999 + 100000*mi + 1000*k)
        bestMdl = train_with_info(modelName, XTrainN, YTrainN, d, sharedInfo);

        if ~sharedInfoLogged
            fprintf("Selected hyperparameters (shared across folds):\n")
            disp(sharedInfo)
            fprintf("Model structure (fold1):\n")
            print_model_structure(bestMdl, modelName)
            sharedInfoLogged = true;
        end

        if modelName == "random_forest"
            yhatTrainN = predict(bestMdl, XTrainN); yhatTrainN = double(yhatTrainN);
            yhatValN   = predict(bestMdl, XValN);   yhatValN   = double(yhatValN);
            yhatTestN  = predict(bestMdl, XTestN);  yhatTestN  = double(yhatTestN);
        else
            yhatTrainN = predict(bestMdl, XTrainN);
            yhatValN   = predict(bestMdl, XValN);
            yhatTestN  = predict(bestMdl, XTestN);
        end

        yhatTrain = double(yhatTrainN) .* sigY + muY;
        yhatVal   = double(yhatValN)   .* sigY + muY;
        yhatTest  = double(yhatTestN)  .* sigY + muY;

        [R2_train(k), MAE_train(k), RMSE_train(k)] = regression_metrics(YTrain, yhatTrain);
        [R2_val(k),   MAE_val(k),   RMSE_val(k)]   = regression_metrics(YVal,   yhatVal);
        [R2_test(k),  MAE_test(k),  RMSE_test(k)]  = regression_metrics(YTest,  yhatTest);

        fprintf("Fold %d metrics: Train R2 %.4f MAE %.4f RMSE %.4f | Val R2 %.4f MAE %.4f RMSE %.4f | Test R2 %.4f MAE %.4f RMSE %.4f\n", ...
            k, R2_train(k), MAE_train(k), RMSE_train(k), R2_val(k), MAE_val(k), RMSE_val(k), R2_test(k), MAE_test(k), RMSE_test(k))

        imp = compute_model_importance(bestMdl, modelName);
        imp_all(k,:) = imp;

        impDat = fullfile(baseDir, sprintf("%s_fold%02d_importance.dat", modelName, k));
        save_importance_dat(impDat, featNames, imp);

        impPng = fullfile(baseDir, sprintf("%s_fold%02d_importance.png", modelName, k));
        plot_importance_bar(impPng, featNames, imp, sprintf("%s | Fold %d | Importance", modelName, k));

        YTrain_all{k} = YTrain;
        YHatTrain_all{k} = yhatTrain;
        YVal_all{k} = YVal;
        YHatVal_all{k} = yhatVal;
        YTest_all{k} = YTest;
        YHatTest_all{k} = yhatTest;

        foldModels{k} = bestMdl;
        foldInfos{k} = sharedInfo;

        mdl = bestMdl;
        info = sharedInfo;

        foldMat = fullfile(baseDir, sprintf("%s_fold%02d_pack.mat", modelName, k));
        save(foldMat, ...
            "mdl", "info", "muX", "sigX", "muY", "sigY", "idxTrain", "idxVal", "idxTest", ...
            "R2_train", "R2_val", "R2_test", "MAE_train", "MAE_val", "MAE_test", "RMSE_train", "RMSE_val", "RMSE_test", ...
            "-v7.3")

        if isvalid(hwb)
            waitbar(min(1, k/totalSteps), hwb, sprintf("Fold %d/%d done | %s", k, K, modelName));
            drawnow limitrate
        end
    end

    save_dat_kfold(fullfile(baseDir, sprintf("%s_train.dat", modelName)), YTrain_all, YHatTrain_all);
    save_dat_kfold(fullfile(baseDir, sprintf("%s_val.dat", modelName)), YVal_all, YHatVal_all);
    save_dat_kfold(fullfile(baseDir, sprintf("%s_test.dat", modelName)), YTest_all, YHatTest_all);

    limParity = compute_parity_lim_all(YTrain_all, YHatTrain_all, YVal_all, YHatVal_all, YTest_all, YHatTest_all);

    plot_parity_kfold(fullfile(baseDir, sprintf("%s_parity_train.png", modelName)), YTrain_all, YHatTrain_all, limParity, modelName+" Train");
    plot_parity_kfold(fullfile(baseDir, sprintf("%s_parity_val.png", modelName)), YVal_all, YHatVal_all, limParity, modelName+" Val");
    plot_parity_kfold(fullfile(baseDir, sprintf("%s_parity_test.png", modelName)), YTest_all, YHatTest_all, limParity, modelName+" Test");

    impMeanK = mean(imp_all, 1);

    impDatK = fullfile(baseDir, sprintf("%s_%dfold_importance.dat", modelName, K));
    save_importance_dat(impDatK, featNames, impMeanK);

    impPngK = fullfile(baseDir, sprintf("%s_%dfold_importance.png", modelName, K));
    plot_importance_bar(impPngK, featNames, impMeanK, sprintf("%s | %d-fold mean | Importance", modelName, K));

    safeCloseWaitbar(hwb);
    clear cleanupObj

    fprintf("\n===== %s | %d-fold average =====\n", modelName, K)
    fprintf("Train: R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_train), mean(MAE_train), mean(RMSE_train))
    fprintf("Val  : R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_val), mean(MAE_val), mean(RMSE_val))
    fprintf("Test : R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_test), mean(MAE_test), mean(RMSE_test))
end

tTotal = toc(tStart);
fprintf("\nTotal runtime: %.2f seconds (%.2f minutes)\n", tTotal, tTotal/60)

diary off

function [mdl, info, yhatValN, yhatTrSubN] = train_one_trial(modelName, XTrainN, YTrainN, XValN, YValN, XTrainSub, YTrainSub, d, mi, k, tr)
rng(42 + 1000*mi + 10*k + tr)

if modelName == "random_forest"
    cfg.NumTrees = randi([400, 1600]);
    cfg.MinLeafSize = randi([1, 22]);
    cfg.NumPredictorsToSample = randi([max(1, floor(sqrt(d))), d]);
    cfg.InBagFraction = 0.60 + 0.38*rand();
    cfg.Seed = 42 + 1000 * mi + 10 * k + tr;

    mdl = TreeBagger( ...
        cfg.NumTrees, XTrainN, YTrainN, ...
        "Method", "regression", ...
        "MinLeafSize", cfg.MinLeafSize, ...
        "NumPredictorsToSample", cfg.NumPredictorsToSample, ...
        "OOBPrediction", "on", ...
        "OOBPredictorImportance", "on", ...
        "InBagFraction", cfg.InBagFraction ...
    );

    yhatValN = predict(mdl, XValN); yhatValN = double(yhatValN);
    yhatTrSubN = predict(mdl, XTrainSub); yhatTrSubN = double(yhatTrSubN);

    info = cfg;
    return
end

if modelName == "xgboost_style"
    cfg.LearnRate = 10.^(-2.2 + (log10(0.12)-log10(0.006))*rand());
    cfg.NumLearningCycles = randi([700, 3600]);
    cfg.MaxDepth = randi([3, 9]);
    cfg.MinLeafSize = randi([1, 28]);
    cfg.Subsample = 0.55 + 0.40*rand();
    cfg.Colsample = 0.55 + 0.40*rand();
    cfg.Seed = 42 + 1000 * mi + 10 * k + tr;

    maxNumSplits = max(3, 2^cfg.MaxDepth - 1);
    numVarsToSample = max(1, min(d, round(cfg.Colsample*d)));

    t = templateTree( ...
        "MaxNumSplits", maxNumSplits, ...
        "MinLeafSize", cfg.MinLeafSize, ...
        "NumVariablesToSample", numVarsToSample);

    mdl = fitrensemble( ...
        XTrainN, YTrainN, ...
        "Method", "LSBoost", ...
        "Learners", t, ...
        "NumLearningCycles", cfg.NumLearningCycles, ...
        "LearnRate", cfg.LearnRate, ...
        "FResample", cfg.Subsample);

    yhatValN = predict(mdl, XValN);
    yhatTrSubN = predict(mdl, XTrainSub);

    info = cfg;
    return
end

if modelName == "lightgbm_style"
    cfg.LearnRate = 10.^(-2.4 + (log10(0.10)-log10(0.004))*rand());
    cfg.NumLearningCycles = randi([1200, 6200]);
    cfg.NumLeaves = randi([15, 255]);
    cfg.MinLeafSize = randi([1, 34]);
    cfg.Subsample = 0.55 + 0.40*rand();
    cfg.Colsample = 0.55 + 0.40*rand();
    cfg.Seed = 42 + 1000 * mi + 10 * k + tr;

    maxDepthApprox = max(2, ceil(log2(double(cfg.NumLeaves)+1)));
    maxNumSplits = max(3, 2^maxDepthApprox - 1);
    numVarsToSample = max(1, min(d, round(cfg.Colsample*d)));

    t = templateTree( ...
        "MaxNumSplits", maxNumSplits, ...
        "MinLeafSize", cfg.MinLeafSize, ...
        "NumVariablesToSample", numVarsToSample);

    mdl = fitrensemble( ...
        XTrainN, YTrainN, ...
        "Method", "LSBoost", ...
        "Learners", t, ...
        "NumLearningCycles", cfg.NumLearningCycles, ...
        "LearnRate", cfg.LearnRate, ...
        "FResample", cfg.Subsample);

    yhatValN = predict(mdl, XValN);
    yhatTrSubN = predict(mdl, XTrainSub);

    info = cfg;
    return
end

error("Unknown model: %s", modelName)
end

function mdl = train_with_info(modelName, XTrainN, YTrainN, d, info)
if modelName == "random_forest"
    mdl = TreeBagger( ...
        info.NumTrees, XTrainN, YTrainN, ...
        "Method", "regression", ...
        "MinLeafSize", info.MinLeafSize, ...
        "NumPredictorsToSample", info.NumPredictorsToSample, ...
        "OOBPrediction", "on", ...
        "OOBPredictorImportance", "on", ...
        "InBagFraction", info.InBagFraction ...
    );
    return
end

if modelName == "xgboost_style"
    maxNumSplits = max(3, 2^info.MaxDepth - 1);
    numVarsToSample = max(1, min(d, round(info.Colsample*d)));

    t = templateTree( ...
        "MaxNumSplits", maxNumSplits, ...
        "MinLeafSize", info.MinLeafSize, ...
        "NumVariablesToSample", numVarsToSample);

    mdl = fitrensemble( ...
        XTrainN, YTrainN, ...
        "Method", "LSBoost", ...
        "Learners", t, ...
        "NumLearningCycles", info.NumLearningCycles, ...
        "LearnRate", info.LearnRate, ...
        "FResample", info.Subsample);
    return
end

if modelName == "lightgbm_style"
    maxDepthApprox = max(2, ceil(log2(double(info.NumLeaves)+1)));
    maxNumSplits = max(3, 2^maxDepthApprox - 1);
    numVarsToSample = max(1, min(d, round(info.Colsample*d)));

    t = templateTree( ...
        "MaxNumSplits", maxNumSplits, ...
        "MinLeafSize", info.MinLeafSize, ...
        "NumVariablesToSample", numVarsToSample);

    mdl = fitrensemble( ...
        XTrainN, YTrainN, ...
        "Method", "LSBoost", ...
        "Learners", t, ...
        "NumLearningCycles", info.NumLearningCycles, ...
        "LearnRate", info.LearnRate, ...
        "FResample", info.Subsample);
    return
end

error("Unknown model: %s", modelName)
end

function print_model_structure(mdl, modelName)
if modelName == "random_forest"
    disp(mdl)
    if isprop(mdl, "Trees")
        fprintf("RandomForest: NumTrees=%d\n", numel(mdl.Trees))
    end
    return
end
disp(mdl)
if isprop(mdl, "ModelParameters")
    disp(mdl.ModelParameters)
end
if isprop(mdl, "Trained")
    try
        L = mdl.Trained;
        if iscell(L) && ~isempty(L)
            disp(L{1})
        end
    catch
    end
end
end

function [R2, RMSE] = r2_rmse(y, yhat)
y = y(:);
yhat = yhat(:);
RMSE = sqrt(mean((y-yhat).^2));
den = sum((y - mean(y)).^2);
if den == 0
    R2 = -inf;
else
    R2 = 1 - sum((y - yhat).^2)/den;
end
end

function [R2, MAE, RMSE] = regression_metrics(y, yhat)
y = y(:);
yhat = yhat(:);
MAE = mean(abs(y-yhat));
RMSE = sqrt(mean((y - yhat).^2));
den = sum((y - mean(y)).^2);
if den == 0
    R2 = NaN;
else
    R2 = 1 - sum((y - yhat).^2) / den;
end
end

function save_dat_kfold(fname, yCells, yhatCells)
K = numel(yCells);
lens = zeros(K, 1);
for k = 1:K
    lens(k) = numel(yCells{k});
end
Lmax = max(lens);
M = NaN(Lmax, 2*K);
for k = 1:K
    y = yCells{k}(:);
    yh = yhatCells{k}(:);
    M(1:numel(y), 2*k-1) = y;
    M(1:numel(yh), 2*k) = yh;
end
fid = fopen(fname, "w");
hdr = '#';
for k = 1:K
    hdr = [hdr, sprintf(' Fold%d_True Fold%d_Pred', k, k)];
end
fprintf(fid, "%s\n", hdr);
fmt = [repmat(' %.8f', 1, 2*K), '\n'];
for i = 1:Lmax
    fprintf(fid, fmt, M(i, :));
end
fclose(fid);
end

function lim = compute_parity_lim_all(YTr, YHTr, YVa, YHVa, YTe, YHTe)
allv = [];
K = numel(YTr);
for k = 1:K
    allv = [allv; YTr{k}(:); YHTr{k}(:); YVa{k}(:); YHVa{k}(:); YTe{k}(:); YHTe{k}(:)];
end
limL = min(allv);
limU = max(allv);
pad = 0.02 * (limU - limL + eps);
lim = [limL - pad, limU + pad];
end

function plot_parity_kfold(saveName, yCells, yhatCells, lim, titleStr)
K = numel(yCells);
fig = figure("Visible", "off");
maxPtsPerFold = inf;
cols = lines(K);
hold on
for k = 1:K
    y = yCells{k}(:);
    yh = yhatCells{k}(:);
    n = numel(y);
    if isfinite(maxPtsPerFold) && n > maxPtsPerFold
        idx = randperm(n, maxPtsPerFold);
        y = y(idx);
        yh = yh(idx);
    end
    plot(y, yh, '.', "Color", cols(k, :), "MarkerSize", 8);
end
plot(lim, lim, "k--", "LineWidth", 1.2);
xlabel("True");
ylabel("Predicted");
title(titleStr);
xlim(lim);
ylim(lim);
axis equal
axis square
grid on
box on
lg = strings(K, 1);
for k = 1:K
    lg(k) = "Fold " + string(k);
end
legend(lg, "Location", "bestoutside");
set(fig, "PaperPositionMode", "auto");
print(fig, saveName, "-dpng", "-r600");
close(fig);
end

function safeCloseWaitbar(h)
if ~isempty(h) && isvalid(h)
    close(h)
end
end

function imp = compute_model_importance(mdl, modelName)
if modelName == "random_forest"
    if isprop(mdl, "OOBPermutedPredictorDeltaError")
        imp = double(mdl.OOBPermutedPredictorDeltaError);
    else
        imp = NaN(1, size(mdl.X,2));
    end
    imp = imp(:).';
    return
end
try
    imp = predictorImportance(mdl);
catch
    imp = NaN(1, numel(mdl.PredictorNames));
end
imp = double(imp(:).');
end

function save_importance_dat(fname, featNames, imp)
[impS, ord] = sort(imp, "descend");
fid = fopen(fname, "w");
fprintf(fid, "Rank Feature Importance\n");
for r = 1:numel(ord)
    j = ord(r);
    fprintf(fid, "%d %s %.10e\n", r, char(featNames(j)), impS(r));
end
fclose(fid);
end

function plot_importance_bar(saveName, featNames, imp, titleStr)
[impS, ord] = sort(imp, "ascend");
fig = figure("Visible", "off");
barh(impS);
yticks(1:numel(ord));
yticklabels(featNames(ord));
xlabel("Importance");
title(titleStr);
grid on
box on
set(fig, "PaperPositionMode", "auto");
print(fig, saveName, "-dpng", "-r600");
close(fig);
end
