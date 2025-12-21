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

T = readtable(dataFile, VariableNamingRule = "preserve");

X = T{:, 6:13};
Y = T{:, end};
X = double(X);
Y = double(Y);

N = size(X, 1);
K = 5; % 5

cv = cvpartition(N, "KFold", K);

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

netCfg.inDim = 8;
netCfg.width = 128; % 32~256, 32 64 ..., high
netCfg.depth = 3; % 1~6, 1 2 ..., high
netCfg.drop = 0.05; % 0~0.3, low
netCfg.glu = true; % true/false

trainCfg.maxEpoch = 200; % 50~400, high
trainCfg.batch = 64; % 16~256, 16 32 ..., high
trainCfg.lr = 1e-3; % 1e-4 ~ 3e-3, high
trainCfg.wd = 3e-4; % 0 ~ 1e-3, high
trainCfg.gradClip = 1.0
trainCfg.lrDropEvery = 40; % 10 ~ 50, high
trainCfg.lrDropFactor = 0.5
trainCfg.patience = 40; % 10 ~ 60, high
trainCfg.emaDecay = 0.995; % 0.95 ~ 0.999, high
trainCfg.ensN = 3; % 1 ~ 10, high

valRatio = 0.2;

totalSteps = K * trainCfg.maxEpoch;
stepNow = 0;

hwb = waitbar(0, "Starting...", "Name", "progress");
cleanupObj = onCleanup(@()safeCloseWaitbar(hwb));

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

    XTrainN = single(XTrainN);
    XValN = single(XValN);
    XTestN = single(XTestN);
    YTrainN = single(YTrainN);
    YValN = single(YValN);
    YTestN = single(YTestN);

    fprintf("Config: width=%d depth=%d drop=%.3f glu=%d lr=%.2e wd=%.2e batch=%d ema=%.4f ens=%d epochs=%d\n", ...
        netCfg.width, netCfg.depth, netCfg.drop, netCfg.glu, trainCfg.lr, trainCfg.wd, trainCfg.batch, trainCfg.emaDecay, trainCfg.ensN, trainCfg.maxEpoch)

    YPredTrainN = zeros(numel(YTrainN), trainCfg.ensN, 'single');
    YPredValN = zeros(numel(YValN), trainCfg.ensN, 'single');
    YPredTestN = zeros(numel(YTestN), trainCfg.ensN, 'single');

    nets = cell(trainCfg.ensN, 1);
    emaNets = cell(trainCfg.ensN, 1);
    bestInfos = cell(trainCfg.ensN, 1);

    for e = 1:trainCfg.ensN
        rng(42+1000*k+e)
        dlnet = init_resglu_dlnet(netCfg);

        if k == 1 && e == 1
            fprintf("\n==== Fold %d/%d | member %d/%d | Network dump (only once) ====\n", k, K, e, trainCfg.ensN);
            log_dlnet_to_log(dlnet);
        end

        [emaNet, bestNet, bestInfo] = train_dlnet_reg(dlnet, XTrainN, YTrainN, XValN, YValN, trainCfg, hwb, k, K, e, trainCfg.ensN, stepNow, totalSteps);
        stepNow = stepNow + bestInfo.epochsRan;

        nets{e} = bestNet;
        emaNets{e} = emaNet;
        bestInfos{e} = bestInfo;

        YPredTrainN(:, e) = predict_dlnet(emaNet, XTrainN, trainCfg.batch);
        YPredValN(:, e) = predict_dlnet(emaNet, XValN, trainCfg.batch);
        YPredTestN(:, e) = predict_dlnet(emaNet, XTestN, trainCfg.batch);

        fprintf("Fold %d ensN %d bestEpoch=%d bestValLoss=%.6f\n", k, e, bestInfo.bestEpoch, bestInfo.bestValLoss)
    end

    yhatTrain = double(mean(YPredTrainN, 2)) .* sigY + muY;
    yhatVal = double(mean(YPredValN, 2)) .* sigY + muY;
    yhatTest = double(mean(YPredTestN, 2)) .* sigY + muY;

    [R2_train(k), MAE_train(k), RMSE_train(k)] = regression_metrics(YTrain, yhatTrain);
    [R2_val(k), MAE_val(k), RMSE_val(k)] = regression_metrics(YVal, yhatVal);
    [R2_test(k), MAE_test(k), RMSE_test(k)] = regression_metrics(YTest, yhatTest);

    fprintf("Fold %d metrics: Train R2 %.4f MAE %.4f RMSE %.4f | Val R2 %.4f MAE %.4f RMSE %.4f | Test R2 %.4f MAE %.4f RMSE %.4f\n", ...
        k, R2_train(k), MAE_train(k), RMSE_train(k), R2_val(k), MAE_val(k), RMSE_val(k), R2_test(k), MAE_test(k), RMSE_test(k))

    YTrain_all{k} = YTrain;
    YHatTrain_all{k} = yhatTrain;
    YVal_all{k} = YVal;
    YHatVal_all{k} = yhatVal;
    YTest_all{k} = YTest;
    YHatTest_all{k} = yhatTest;

    foldMat = fullfile(baseDir, sprintf("fold%02d_pack.mat", k));
    save(foldMat, ...
        "nets", "emaNets", "bestInfos", "muX", "sigX", "muY", "sigY", "idxTrain", "idxVal", "idxTest", ...
        "netCfg", "trainCfg", "R2_train", "R2_val", "R2_test", "MAE_train", "MAE_val", "MAE_test", "RMSE_train", "RMSE_val", "RMSE_test", ...
        "-v7.3")
end

save_dat_kfold(fullfile(baseDir, "train.dat"), YTrain_all, YHatTrain_all);
save_dat_kfold(fullfile(baseDir, "val.dat"), YVal_all, YHatVal_all);
save_dat_kfold(fullfile(baseDir, "test.dat"), YTest_all, YHatTest_all);

[limParity] = compute_parity_lim_all(YTrain_all, YHatTrain_all, YVal_all, YHatVal_all, YTest_all, YHatTest_all);

plot_parity_kfold(fullfile(baseDir, "parity_train.png"), YTrain_all, YHatTrain_all, limParity, "Train");
plot_parity_kfold(fullfile(baseDir, "parity_val.png"), YVal_all, YHatVal_all, limParity, "Val");
plot_parity_kfold(fullfile(baseDir, "parity_test.png"), YTest_all, YHatTest_all, limParity, "Test");

safeCloseWaitbar(hwb);
clear cleanupObj

fprintf("\n===== %d-fold average =====\n", K)
fprintf("Train: R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_train), mean(MAE_train), mean(RMSE_train))
fprintf("Val  : R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_val), mean(MAE_val), mean(RMSE_val))
fprintf("Test : R2 = %.4f, MAE = %.4f, RMSE = %.4f\n", mean(R2_test), mean(MAE_test), mean(RMSE_test))

tTotal = toc(tStart);
fprintf("\nTotal runtime: %.2f seconds (%.2f minutes)\n", tTotal, tTotal/60)

diary off

function dlnet = init_resglu_dlnet(cfg)
in = featureInputLayer(cfg.inDim, Normalization = "none", Name = "in");

stem = [; ...
    fullyConnectedLayer(cfg.width, Name = "stem_fc"); ...
    layerNormalizationLayer(Name = "stem_ln"); ...
    swishLayer(Name = "stem_sw"); ...
    dropoutLayer(cfg.drop, Name = "stem_dp"); ...
    ];

lgraph = layerGraph([in; stem]);

for i = 1:cfg.depth
    tag = "b" + string(i);

    if cfg.glu
        b = [; ...
            fullyConnectedLayer(cfg.width*2, Name = tag+"_fc1"); ...
            layerNormalizationLayer(Name = tag+"_ln1"); ...
            swishLayer(Name = tag+"_sw1"); ...
            dropoutLayer(cfg.drop, Name = tag+"_dp1"); ...
            functionLayer(@(X)glu_split(X), Formattable = true, Name = tag+"_glu"); ...
            fullyConnectedLayer(cfg.width, Name = tag+"_fc2"); ...
            layerNormalizationLayer(Name = tag+"_ln2"); ...
            additionLayer(2, Name = tag+"_add"); ...
            swishLayer(Name = tag+"_sw2"); ...
            dropoutLayer(cfg.drop, Name = tag+"_dp2"); ...
            ];
    else
        b = [; ...
            fullyConnectedLayer(cfg.width, Name = tag+"_fc1"); ...
            layerNormalizationLayer(Name = tag+"_ln1"); ...
            swishLayer(Name = tag+"_sw1"); ...
            dropoutLayer(cfg.drop, Name = tag+"_dp1"); ...
            fullyConnectedLayer(cfg.width, Name = tag+"_fc2"); ...
            layerNormalizationLayer(Name = tag+"_ln2"); ...
            additionLayer(2, Name = tag+"_add"); ...
            swishLayer(Name = tag+"_sw2"); ...
            dropoutLayer(cfg.drop, Name = tag+"_dp2"); ...
            ];
    end

    lgraph = addLayers(lgraph, b);

    if i == 1
        inName = "stem_dp";
    else
        inName = "b" + string(i-1) + "_dp2";
    end

    lgraph = connectLayers(lgraph, inName, tag+"_fc1");
    lgraph = connectLayers(lgraph, inName, tag+"_add/in2");
end

tail = [; ...
    fullyConnectedLayer(64, Name = "tail_fc1"); ...
    layerNormalizationLayer(Name = "tail_ln1"); ...
    swishLayer(Name = "tail_sw1"); ...
    dropoutLayer(max(0.02, cfg.drop/2), Name = "tail_dp1"); ...
    fullyConnectedLayer(32, Name = "tail_fc2"); ...
    swishLayer(Name = "tail_sw2"); ...
    fullyConnectedLayer(1, Name = "out_fc"); ...
    ];

lgraph = addLayers(lgraph, tail);

if cfg.depth >= 1
    last = "b" + string(cfg.depth) + "_dp2";
else
    last = "stem_dp";
end

lgraph = connectLayers(lgraph, last, "tail_fc1");

dlnet = dlnetwork(lgraph);
end

function Y = glu_split(X)
C = size(X, 1);
h = C / 2;
A = X(1:h, :, :, :);
B = X(h+1:end, :, :, :);
Y = A .* sigmoid(B);
end

function [emaNet, bestNet, info] = train_dlnet_reg(dlnet, XTrain, YTrain, XVal, YVal, cfg, hwb, k, K, e, ensN, stepBase, totalSteps)
mb = cfg.batch;
numObs = size(XTrain, 1);
numIter = ceil(numObs/mb);

trailingAvg = [];
trailingAvgSq = [];

bestVal = inf;
bestNet = dlnet;
bestEpoch = 0;

emaNet = dlnet;
emaDecay = cfg.emaDecay;

pat = 0;
epochsRan = 0;

for epoch = 1:cfg.maxEpoch
    epochsRan = epochsRan + 1;

    if isvalid(hwb)
        prog = min(1, (stepBase + epoch)/max(1, totalSteps));
        msg = sprintf("Fold %d/%d | ensN %d/%d | Epoch %d/%d", k, K, e, ensN, epoch, cfg.maxEpoch);
        waitbar(prog, hwb, msg);
        drawnow limitrate
    end

    lr = cfg.lr * (cfg.lrDropFactor^floor((epoch - 1)/cfg.lrDropEvery));

    idx = randperm(numObs);

    for it = 1:numIter
        i1 = (it - 1) * mb + 1;
        i2 = min(it*mb, numObs);
        bi = idx(i1:i2);

        Xb = dlarray(XTrain(bi, :)', "CB");
        Yb = dlarray(YTrain(bi, :)', "CB");

        [loss, grads] = dlfeval(@modelGradients, dlnet, Xb, Yb, cfg.wd);

        grads = dlupdate(@(g)clipGrad(g, cfg.gradClip), grads);

        [dlnet, trailingAvg, trailingAvgSq] = adamupdate(dlnet, grads, trailingAvg, trailingAvgSq, epoch*100000+it, lr);

        emaNet = emaUpdate(emaNet, dlnet, emaDecay);
    end

    valPred = predict_dlnet(emaNet, XVal, cfg.batch);
    vLoss = mean((valPred - YVal).^2, "all");

    if vLoss < bestVal
        bestVal = vLoss;
        bestNet = dlnet;
        bestEpoch = epoch;
        pat = 0;
    else
        pat = pat + 1;
    end

    if pat >= cfg.patience
        break
    end
end

info.bestValLoss = bestVal;
info.bestEpoch = bestEpoch;
info.epochsRan = epochsRan;
end

function [loss, grads] = modelGradients(dlnet, X, Y, wd)
YHat = forward(dlnet, X);
mse = mean((YHat - Y).^2, "all");

l2 = 0;
L = dlnet.Learnables;

for i = 1:size(L, 1)
    p = L.Value{i};
    if isnumeric(p) || isa(p, "dlarray")
        l2 = l2 + sum(p.^2, "all");
    end
end

loss = mse + wd * l2;

grads = dlgradient(loss, dlnet.Learnables);
end

function g = clipGrad(g, thr)
if isempty(g)
    return
end
n = sqrt(sum(g.^2, "all"));
if n > thr
    g = g * (thr / (n + 1e-12));
end
end

function emaNet = emaUpdate(emaNet, net, decay)
L = emaNet.Learnables;
for i = 1:size(L, 1)
    v1 = L.Value{i};
    v2 = net.Learnables.Value{i};
    L.Value{i} = decay * v1 + (1 - decay) * v2;
end
emaNet.Learnables = L;
end

function yhat = predict_dlnet(net, X, batch)
n = size(X, 1);
yhat = zeros(n, 1, 'single');
i = 1;
while i <= n
    j = min(i+batch-1, n);
    Xb = dlarray(single(X(i:j, :))', "CB");
    yb = extractdata(forward(net, Xb))';
    yhat(i:j) = single(yb);
    i = j + 1;
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

function log_dlnet_to_log(dlnet)
layers = dlnet.Layers;
fprintf("---- Layers: %d ----\n", numel(layers));
for i = 1:numel(layers)
    lyr = layers(i);
    fprintf("%3d | %-28s | %s\n", i, string(lyr.Name), class(lyr));
end

L = dlnet.Learnables;
fprintf("---- Learnables: %d ----\n", size(L, 1));

totalParams = 0;
for i = 1:size(L, 1)
    v = L.Value{i};

    sz = size(v);
    n = numel(v);

    totalParams = totalParams + n;

    fprintf("%3d | %-24s | %-18s | size=%s | numel=%d\n", ...
        i, string(L.Layer(i)), string(L.Parameter(i)), mat2str(sz), n);
end

fprintf("Total learnable parameters: %d\n", totalParams);
fprintf("==== End network dump ====\n\n");
end
