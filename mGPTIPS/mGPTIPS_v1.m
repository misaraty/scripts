clear;clc;close all;
scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir);scriptDir=pwd;end
cd(scriptDir);
SEED=42;
rng(SEED,"twister");
widColon="MATLAB:colon:operandsMustBeRealScalars";
wColon=warning("query",widColon);
warning("off",widColon);
XLSX="data_v4.xlsx";
OUTDIR="gptips2_bandgap_split_out";
if ~exist(OUTDIR,"dir");mkdir(OUTDIR);end
TEST_RATIO=0.2;
VAL_RATIO=0.2;  % 0.2*(1-TEST_RATIO)
GP.pop_size=300;
GP.num_gen=150;
GP.max_depth=5;
GP.max_genes=8;
GP.funcs={'times','minus','plus','rdivide'};
% GP.funcs={'times','minus','plus','rdivide','square','sin','cos','exp','mult3','add3','sqrt','cube','power','negexp','neg','abs','log'};
GP.const_range=[-5 5];
GP.num_runs=20;
GP.verbose=10;
GP.elite_fraction=0.2;
GP.tournament_size=10;
GP.use_lexicographic=true;
GP.cache_fitness=true;
GP.use_zscore_x=true;
GP.use_zscore_y=true;
GP.y_outlier_clip_sigma=6;
GP.parsimony_on=true;
GP.parsimony_p_pareto=0.35;
GP.timeout_sec=0;
GP.select_by="pareto_valrmse_complexity";
GP.parity_dpi=600;
GP.marker_size=36;
GP.line_width=1.2;
GP.fig_w_in=6.0;
GP.fig_h_in=6.0;
featNames=["r_A","r_B'","r_B''","r_X","χ_A","χ_B'","χ_B''","χ_X"];
% featNames=["r_A","r_B'","r_B''","r_X"];
% featNames=["χ_A","χ_B'","χ_B''","χ_X"];
targetName="Bandgap";
T=readtable(XLSX,"VariableNamingRule","preserve");
[X,y,usedFeatNames]=extractXY_from_table(T,featNames,targetName);
mask=all(isfinite(X),2)&isfinite(y);
X=X(mask,:);y=y(mask,:);
N=size(X,1);
fprintf("Loaded %d samples, %d features.\n",N,size(X,2));
disp("Features used:");disp(usedFeatNames);
idx=randperm(N);
nTest=max(1,round(TEST_RATIO*N));
idxTest=idx(1:nTest);
idxRest=idx(nTest+1:end);
X_te=X(idxTest,:);y_te=y(idxTest,:);
X_rest=X(idxRest,:);y_rest=y(idxRest,:);
nRest=size(X_rest,1);
idx2=randperm(nRest);
nVal=max(1,round(VAL_RATIO*nRest));
idxVal_local=idx2(1:nVal);
idxTr_local=idx2(nVal+1:end);
X_tr=X_rest(idxTr_local,:);y_tr=y_rest(idxTr_local,:);
X_va=X_rest(idxVal_local,:);y_va=y_rest(idxVal_local,:);
[X_tr_z,muX,sigX]=zscore_safe(X_tr,GP.use_zscore_x);
X_va_z=(X_va-muX)./sigX;
X_te_z=(X_te-muX)./sigX;
[y_tr_z,muY,sigY]=zscore_safe(y_tr,GP.use_zscore_y);
y_va_z=(y_va-muY)./sigY;
y_te_z=(y_te-muY)./sigY;
y_tr_z=clip_sigma(y_tr_z,GP.y_outlier_clip_sigma);
gp=run_gptips2_multigene_regression(X_tr_z,y_tr_z,X_va_z,y_va_z,X_te_z,y_te_z,usedFeatNames,GP,SEED);
candTags=["valbest","best"];
chosenTag="";
if isfield(gp,"results")
    for i=1:numel(candTags)
        if has_model_tag(gp,candTags(i))
            chosenTag=string(candTags(i));
            break;
        end
    end
end
if strlength(chosenTag)==0
    chosenTag="best";
end
bestParetoTag=pick_model_by_pareto(gp,X_tr_z,y_tr_z,X_va_z,y_va_z,candTags,GP);
if strlength(bestParetoTag)>0
    modelTag=bestParetoTag;
else
    modelTag=chosenTag;
end
fprintf("Model selected: %s\n",modelTag);
[pred_tr_z,pred_va_z,pred_te_z]=predict_with_gp(gp,modelTag,X_tr_z,X_va_z,X_te_z);
pred_tr=pred_tr_z*sigY+muY;
pred_va=pred_va_z*sigY+muY;
pred_te=pred_te_z*sigY+muY;
m_tr=calc_metrics(y_tr,pred_tr);
m_va=calc_metrics(y_va,pred_va);
m_te=calc_metrics(y_te,pred_te);
fprintf("Train: R2=%.4f  MAE=%.4f  RMSE=%.4f\n",m_tr(1),m_tr(2),m_tr(3));
fprintf("Val  : R2=%.4f  MAE=%.4f  RMSE=%.4f\n",m_va(1),m_va(2),m_va(3));
fprintf("Test : R2=%.4f  MAE=%.4f  RMSE=%.4f\n",m_te(1),m_te(2),m_te(3));
rawEqn=extract_equation_raw_silent(gp,modelTag);
complexity=extract_complexity(gp,modelTag);
fh=extract_fitness_history(gp);
runDat=struct();
runDat.seed=SEED;
runDat.usedFeatNames=usedFeatNames;
runDat.idxTrain=idxRest(idxTr_local);
runDat.idxVal=idxRest(idxVal_local);
runDat.idxTest=idxTest;
runDat.muX=muX;runDat.sigX=sigX;
runDat.muY=muY;runDat.sigY=sigY;
runDat.y_tr=y_tr;runDat.pred_tr=pred_tr;
runDat.y_va=y_va;runDat.pred_va=pred_va;
runDat.y_te=y_te;runDat.pred_te=pred_te;
runDat.metrics_train=m_tr;
runDat.metrics_val=m_va;
runDat.metrics_test=m_te;
runDat.modelTag=modelTag;
runDat.eqn_raw=rawEqn;
runDat.complexity=complexity;
runDat.fitness_history=fh;
save(fullfile(OUTDIR,"split_results.mat"),"-struct","runDat");
write_parity_dat(y_tr,pred_tr,fullfile(OUTDIR,"parity_train.dat"));
write_parity_dat(y_va,pred_va,fullfile(OUTDIR,"parity_val.dat"));
write_parity_dat(y_te,pred_te,fullfile(OUTDIR,"parity_test.dat"));
make_parity_plot(y_tr,pred_tr,"Train Parity",fullfile(OUTDIR,"parity_train.png"),GP,5);
make_parity_plot(y_va,pred_va,"Val Parity",fullfile(OUTDIR,"parity_val.png"),GP,1);
make_parity_plot(y_te,pred_te,"Test Parity",fullfile(OUTDIR,"parity_test.png"),GP,9);
if ~isempty(fh)
    writematrix([(1:numel(fh))',fh(:)],fullfile(OUTDIR,"fitness.dat"),"Delimiter","\t");
    make_fitness_plot(fh,"Fitness (lower=better)",fullfile(OUTDIR,"fitness.png"),GP);
end
writecell({"Split","R2","MAE","RMSE"},fullfile(OUTDIR,"metrics.dat"),"Delimiter","\t");
writecell({"Train",m_tr(1),m_tr(2),m_tr(3)},fullfile(OUTDIR,"metrics.dat"),"Delimiter","\t","WriteMode","append");
writecell({"Val",m_va(1),m_va(2),m_va(3)},fullfile(OUTDIR,"metrics.dat"),"Delimiter","\t","WriteMode","append");
writecell({"Test",m_te(1),m_te(2),m_te(3)},fullfile(OUTDIR,"metrics.dat"),"Delimiter","\t","WriteMode","append");
fid=fopen(fullfile(OUTDIR,"summary.txt"),"w");
fprintf(fid,"Seed: %d\n\n",SEED);
fprintf(fid,"Train/Val/Test sizes: %d / %d / %d\n\n",numel(runDat.idxTrain),numel(runDat.idxVal),numel(runDat.idxTest));
fprintf(fid,"Selected model: %s\n\n",modelTag);
fprintf(fid,"Train: R2=%.6f  MAE=%.6f  RMSE=%.6f\n",m_tr(1),m_tr(2),m_tr(3));
fprintf(fid,"Val  : R2=%.6f  MAE=%.6f  RMSE=%.6f\n",m_va(1),m_va(2),m_va(3));
fprintf(fid,"Test : R2=%.6f  MAE=%.6f  RMSE=%.6f\n\n",m_te(1),m_te(2),m_te(3));
fprintf(fid,"Complexity: %g\n\n",complexity);
fprintf(fid,"Equation (raw):\n%s\n",rawEqn);
fclose(fid);
fprintf("\nDone. Outputs saved in: %s\n",fullfile(scriptDir,OUTDIR));
warning(wColon.state,widColon);

function [X,y,usedFeatNames]=extractXY_from_table(T,featNames,targetName)
usedFeatNames=featNames;
allVars=string(T.Properties.VariableNames);
for i=1:numel(featNames)
    if ~any(allVars==featNames(i))
        cand=allVars(contains(allVars,erase(featNames(i),"'"),"IgnoreCase",true));
        if isempty(cand)
            error("Feature column not found: %s. Available vars:\n%s",featNames(i),strjoin(allVars,", "));
        else
            usedFeatNames(i)=cand(1);
        end
    end
end
if ~any(allVars==targetName)
    cand=allVars(contains(allVars,targetName,"IgnoreCase",true));
    if isempty(cand)
        error("Target column not found: %s. Available vars:\n%s",targetName,strjoin(allVars,", "));
    else
        targetName=cand(1);
    end
end
n=height(T);
p=numel(usedFeatNames);
X=nan(n,p);
for j=1:p
    col=T.(usedFeatNames(j));
    X(:,j)=col2double(col,usedFeatNames(j));
end
y=col2double(T.(targetName),targetName);
end

function [Z,mu,sig]=zscore_safe(X,doit)
X=double(X);
if ~doit
    mu=zeros(1,size(X,2));
    sig=ones(1,size(X,2));
    Z=X;
    return;
end
if isvector(X)
    mu=mean(X,"omitnan");
    sig=std(X,0,"omitnan");
    if sig==0||~isfinite(sig);sig=1;end
    Z=(X-mu)./sig;
else
    mu=mean(X,1,"omitnan");
    sig=std(X,0,1,"omitnan");
    sig(sig==0|~isfinite(sig))=1;
    Z=(X-mu)./sig;
end
end

function y=clip_sigma(y,nsig)
if isempty(y)||~isfinite(nsig)||nsig<=0;return;end
m=mean(y,"omitnan");
s=std(y,0,"omitnan");
if s==0||~isfinite(s);return;end
lo=m-nsig*s;hi=m+nsig*s;
y(y<lo)=lo;
y(y>hi)=hi;
end

function gp=run_gptips2_multigene_regression(Xtr,ytr,Xva,yva,Xte,yte,featNames,GP,seed)
rng(seed,"twister");
gp=rungp(@configFcn);
    function gp=configFcn(gp)
        gp=gpdefaults();
        gp.userdata.xtrain=Xtr;
        gp.userdata.ytrain=ytr(:);
        gp.userdata.xval=Xva;
        gp.userdata.yval=yva(:);
        gp.userdata.xtest=Xte;
        gp.userdata.ytest=yte(:);
        gp.runcontrol.pop_size=GP.pop_size;
        gp.runcontrol.num_gen=GP.num_gen;
        gp.runcontrol.num_runs=GP.num_runs;
        gp.runcontrol.verbose=GP.verbose;
        if isfield(GP,"timeout_sec") && GP.timeout_sec>0
            try;gp.runcontrol.timeout=GP.timeout_sec;catch;end
        end
        try;gp.runcontrol.parallel.auto=false;catch;end
        try;gp.selection.tournament.size=GP.tournament_size;catch;end
        try;gp.selection.elite_fraction=GP.elite_fraction;catch;end
        try;gp.selection.lexicographic=GP.use_lexicographic;catch;end
        try;gp.fitness.cache=GP.cache_fitness;catch;end
        if isfield(GP,"parsimony_on") && GP.parsimony_on
            try;gp.selection.tournament.p_pareto=GP.parsimony_p_pareto;catch;end
            try;gp.selection.tournament.pareto=GP.parsimony_p_pareto;catch;end
        end
        gp.fitness.fitfun=@regressmulti_fitfun;
        gp.fitness.minimisation=true;
        gp.genes.multigene=true;
        gp.genes.max_genes=GP.max_genes;
        gp.treedef.max_depth=GP.max_depth;
        gp.nodes.inputs.num_inp=size(Xtr,2);
		fn=GP.funcs;
		if isstring(fn);fn=cellstr(fn);end
		if iscell(fn)
			fn=cellfun(@char,fn,"UniformOutput",false);
		else
			fn=cellstr(string(fn));
			fn=cellfun(@char,fn,"UniformOutput",false);
		end
		gp.nodes.functions.name=fn;
        gp.nodes.const.range=GP.const_range;
        try;gp.nodes.terminalSet=cellstr(featNames);catch;end
        try;gp.userdata.user_fcn=@regressmulti_fitfun_validate;catch;end
    end
end

function tf=has_model_tag(gp,tag)
tf=false;
try
    tf=isfield(gp,"results")&&isfield(gp.results,char(tag));
    if tf;return;end
catch
end
try
    runtree(gp,char(tag));
    tf=true;
catch
    tf=false;
end
end

function [ptr,pva,pte]=predict_with_gp(gp,tag,Xtr,Xva,Xte)
tag=char(string(tag));
ptr=[];pva=[];pte=[];
try
    ptr=runtree(gp,tag,Xtr);
    pva=runtree(gp,tag,Xva);
    pte=runtree(gp,tag,Xte);
    return;
catch
end
tmpName=char("tmp_gp_model_"+string(tag)+"_"+string(randi(1e9)));
here=pwd;
tmpM=fullfile(here,[tmpName '.m']);
if isfile(tmpM);delete(tmpM);end
try
    warning('off','MATLAB:colon:operandsMustBeRealScalars');
catch
end
try
    gpmodel2mfile(gp,tag,tmpName);
    addpath(here);
    f=str2func(tmpName);
    ptr=f(Xtr);pva=f(Xva);pte=f(Xte);
    rmpath(here);
    if isfile(tmpM);delete(tmpM);end
catch ME
    if isfile(tmpM);delete(tmpM);end
    error("Prediction failed."+newline+string(ME.message));
end
end

function m=calc_metrics(y,yhat)
y=y(:);yhat=yhat(:);
ssRes=sum((y-yhat).^2);
ssTot=sum((y-mean(y)).^2);
R2=1-ssRes/ssTot;
MAE=mean(abs(y-yhat));
RMSE=sqrt(mean((y-yhat).^2));
m=[R2,MAE,RMSE];
end

function tag=pick_model_by_pareto(gp,Xtr,ytr,Xva,yva,candTags,GP)
tag="";
bestScore=Inf;
bestC=Inf;
for i=1:numel(candTags)
    t=string(candTags(i));
    if ~has_model_tag(gp,t);continue;end
    try
        pva=runtree(gp,char(t),Xva);
        m=calc_metrics(yva,pva);
        rmse=m(3);
        c=extract_complexity(gp,t);
        score=rmse+1e-4*c;
        if score<bestScore || (abs(score-bestScore)<1e-12 && c<bestC)
            bestScore=score;
            bestC=c;
            tag=t;
        end
    catch
    end
end
if strlength(tag)==0
    tag="";
end
end

function make_parity_plot(y,yhat,ttl,savePath,GP,colorIndex)
y=y(:);yhat=yhat(:);
f=figure("Visible","off","Units","inches","Position",[1 1 GP.fig_w_in GP.fig_h_in]);
ax=axes(f);
hold(ax,"on");grid(ax,"on");box(ax,"on");
tc=tab10();
ci=mod(colorIndex-1,size(tc,1))+1;
c=tc(ci,:);
scatter(ax,y,yhat,GP.marker_size,"filled","MarkerFaceColor",c,"MarkerEdgeColor",c,"MarkerFaceAlpha",0.90);
mn=min([y;yhat]);mx=max([y;yhat]);
pad=0.02*(mx-mn);
if ~isfinite(pad)||pad==0;pad=1;end
mn=mn-pad;mx=mx+pad;
plot(ax,[mn mx],[mn mx],"-","LineWidth",GP.line_width,"Color",[0 0 0]);
xlim(ax,[mn mx]);ylim(ax,[mn mx]);
axis(ax,"square");
xlabel(ax,"True");
ylabel(ax,"Predicted");
title(ax,ttl,"Interpreter","none");
set(ax,"FontSize",12,"LineWidth",1);
exportgraphics(f,savePath,"Resolution",GP.parity_dpi);
close(f);
end

function make_fitness_plot(hist,ttl,savePath,GP)
hist=hist(:);
f=figure("Visible","off","Units","inches","Position",[1 1 GP.fig_w_in GP.fig_h_in]);
ax=axes(f);
hold(ax,"on");grid(ax,"on");box(ax,"on");
tc=tab10();
c=tc(2,:);
plot(ax,hist,"-","LineWidth",GP.line_width,"Color",c);
xlabel(ax,"Generation");
ylabel(ax,"Fitness");
title(ax,ttl,"Interpreter","none");
set(ax,"FontSize",12,"LineWidth",1);
exportgraphics(f,savePath,"Resolution",GP.parity_dpi);
close(f);
end

function fh=extract_fitness_history(gp)
fh=[];
candidates={"results.bestfithist","results.stats.bestfitness","results.stats.bestfit","state.bestfitness","state.bestfit","results.fithistory"};
for i=1:numel(candidates)
    v=get_nested_field(gp,candidates{i});
    if isnumeric(v)&&~isempty(v)
        fh=v(:);
        return;
    end
end
end

function v=get_nested_field(s,path)
v=[];
try
    parts=split(string(path),".");
    vtmp=s;
    for k=1:numel(parts)
        p=parts(k);
        if isstruct(vtmp)&&isfield(vtmp,p)
            vtmp=vtmp.(p);
        else
            v=[];
            return;
        end
    end
    v=vtmp;
catch
    v=[];
end
end

function rawEqn=extract_equation_raw_silent(gp,tag)
rawEqn="(Equation extraction via gppretty failed.)";
wid="MATLAB:colon:operandsMustBeRealScalars";
w0=[];
try
    w0=warning("query",wid);
    warning("off",wid);
catch
end
txt="";
try
    txt=evalc("gppretty(gp, char(tag));");
catch
    txt="";
end
try
    if ~isempty(w0);warning(w0.state,wid);end
catch
end
if strlength(string(txt))==0
    rawEqn="(Equation extraction via gppretty failed.)";
else
    s=string(txt);
    s=regexprep(s,"\[[^\]]*\]","");
    s=regexprep(s,"\x08","");
    s=regexprep(s,"\r","");
    lines=splitlines(s);
    keep=true(size(lines));
    for i=1:numel(lines)
        L=strtrim(lines(i));
        if startsWith(L,"Warning:","IgnoreCase",true) || startsWith(L,"警告:") || contains(L,"位置:") || contains(L,"opentoline(") || contains(L,"errorDocCallback(")
            keep(i)=false;
        end
    end
    lines=lines(keep);
    rawEqn=strtrim(strjoin(lines,newline));
    if strlength(rawEqn)==0
        rawEqn="(Equation extraction via gppretty failed.)";
    end
end
end

function c=extract_complexity(gp,tag)
c=NaN;
candidates={"results.(tag).complexity","results.(tag).comp","results.best.complexity"};
for i=1:numel(candidates)
    path=replace(candidates{i},"(tag)",string(tag));
    v=get_nested_field(gp,path);
    if isnumeric(v)&&~isempty(v)
        c=double(v(1));
        return;
    end
end
try
    rawEqn=string(evalc("gppretty(gp, char(tag));"));
    c=strlength(rawEqn);
catch
    c=NaN;
end
end

function v=col2double(col,colName)
try
    if isnumeric(col);v=double(col);return;end
    if islogical(col);v=double(col);return;end
    if iscategorical(col);col=string(col);end
    if isstring(col);v=str2double(col);return;end
    if iscell(col)
        if all(cellfun(@(x) isnumeric(x)&&isscalar(x),col))
            v=cellfun(@double,col);
            return;
        else
            colStr=string(col);
            v=str2double(colStr);
            return;
        end
    end
    if ischar(col);v=str2double(string(col));return;end
    v=double(col);
catch ME
    error("Failed to convert column '%s' to numeric. MATLAB type: %s. Error: %s",colName,class(col),ME.message);
end
end

function write_parity_dat(y,yhat,savePath)
fid=fopen(savePath,"w");
fprintf(fid,"# True\tPredicted");
fclose(fid);
M=[y(:),yhat(:)];
writematrix(M,savePath,"Delimiter","\t","WriteMode","append");
end

function C=tab10()
C=[0.1216 0.4667 0.7059;1.0000 0.4980 0.0549;0.1725 0.6275 0.1725;0.8392 0.1529 0.1569;0.5804 0.4039 0.7412;0.5490 0.3373 0.2941;0.8902 0.4667 0.7608;0.4980 0.4980 0.4980;0.7373 0.7412 0.1333;0.0902 0.7451 0.8118];
end
