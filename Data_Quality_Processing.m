%% ================================================================
%  Influence_and_Robustness.m
%  刀切法影响点诊断 + 偏相关 + 双口径对照
%% ================================================================
clear; clc;

%% ---- 载入（沿用主脚本的按名对齐）----
A = load('Individual_wPLI_Data_TrainRest_Alpha.mat');
B = load('Individual_wPLI_Data_RAT1RAT2_Alpha.mat');
assert(isequal(A.final_subs, B.final_subs));
subs = A.final_subs(:);  cl = A.chan_labels;  n = numel(subs);

CT_subs = {'CT001','CT004','CT016','CT022','CT034','CT035','CT042','CT055','CT061','CT071', ...
           'CT077','CT080','CT090','CT094','CT096','CT103','CT105','CT106','CT113','CT116'};
FA_subs = {'FA002','FA003','FA006','FA009','FA010','FA012','FA013','FA024','FA025','FA033', ...
           'FA050','FA059','FA062','FA063','FA070','FA075','FA078','FA092','FA098','FA101', ...
           'FA112','FA118','FA119'};
dCT = [33,24,25,26,25,25,30,24,23,33,23,19,31,23,24,31,29,28,33,32]' - ...
      [26,23,23,17,25,26,18,23,29,31,22,24,26,19,29,29,28,27,29,25]';
dFA = [22,32,32,30,20,33,24,24,15,35,19,23,24,20,26,32,21,25,22,23,28,25,27]' - ...
      [20,27,29,24,18,31,29,13,13,23, 8,24,20,21,13,27,20,25,23,23,28,22,25]';
BEH  = containers.Map([CT_subs,FA_subs], num2cell([dCT;dFA]));
dRAT = cellfun(@(s) BEH(s), subs);
isFA = startsWith(subs,'FA');

E   = @(D,a,b) squeeze(D(:, strcmp(cl,a), strcmp(cl,b)));
roi = mean([E(A.final_wPLI_DIFF,'Pz','Fp1'), E(A.final_wPLI_DIFF,'Pz','Fp2'), ...
            E(A.final_wPLI_DIFF,'Pz','F3'),  E(A.final_wPLI_DIFF,'Pz','F4')], 2);
f3  = E(A.final_wPLI_DIFF,'Pz','F3');
f7  = E(B.final_wPLI_DIFF,'F7','FC1');

TESTS = { 'H1_ROI',       roi, dRAT;
          'H1_PzF3',      f3,  dRAT;
          'H2_F7FC1',     f7,  dRAT;
          'H3_F3xF7FC1',  f3,  f7   };

%% ================================================================
%  [A] 刀切法：逐一剔除，看谁在撑效应
%% ================================================================
fprintf('\n===== [A] Leave-one-out 影响点诊断（FA 组）=====\n');
flag = {'FA006','FA013','FA070','FA078'};       % >20% 通道剔除者
for t = 1:size(TESTS,1)
    x = TESTS{t,2}(isFA);  y = TESTS{t,3}(isFA);
    s = subs(isFA);  m = numel(s);
    r_full = corr(x,y);
    r_loo = nan(m,1);
    for i = 1:m, k = true(m,1); k(i) = false; r_loo(i) = corr(x(k),y(k)); end
    d = r_loo - r_full;                          % 负=该人在支撑效应
    [~,ord] = sort(abs(d),'descend');
    fprintf('\n  %s : 全样本 r = %+.3f (n=%d)\n', TESTS{t,1}, r_full, m);
    fprintf('    最具影响力的 5 人（剔除后 r 及变化量）:\n');
    for j = 1:5
        i = ord(j);
        mark = ''; if ismember(s{i},flag), mark = '  ← >20%通道'; end
        fprintf('      %-7s r=%+.3f  Δ=%+.3f%s\n', s{i}, r_loo(i), d(i), mark);
    end
    fprintf('    LOO 区间: [%+.3f, %+.3f]', min(r_loo), max(r_loo));
    if sign(min(r_loo))~=sign(max(r_loo)) || min(abs(r_loo))<0.2
        fprintf('   ⚠️ 不稳定');
    end
    fprintf('\n');
    % 同时剔除 4 名高剔除率被试
    k4 = ~ismember(s, flag);
    [r4,p4] = corr(x(k4), y(k4));
    fprintf('    剔除全部 4 名高剔除率被试: r = %+.3f, p = %.3f (n=%d)\n', r4, p4, sum(k4));
end

%% ================================================================
%% [B-fix] 正确的 bootstrap + Fisher CI 对照
rng(2024);
for t = 1:size(TESTS,1)
    x = TESTS{t,2}(isFA);  y = TESTS{t,3}(isFA);  m = numel(x);
    [rp,pp] = corr(x,y);
    [rs,ps] = corr(x,y,'Type','Spearman');

    bs = nan(5000,1);
    for b = 1:5000
        k = randi(m, m, 1);          % ★ 同一套索引
        if std(x(k))>0 && std(y(k))>0, bs(b) = corr(x(k), y(k)); end
    end
    ci_b = prctile(bs(~isnan(bs)), [2.5 97.5]);
    ci_f = tanh(atanh(rp) + [-1 1]*1.96/sqrt(m-3));   % Fisher 参数 CI

    fprintf('  %-14s Pearson %+.3f (p=%.3f) | Spearman %+.3f (p=%.3f)\n', ...
            TESTS{t,1}, rp,pp, rs,ps);
    fprintf('  %-14s boot 95%%CI [%+.3f,%+.3f]%s | Fisher 95%%CI [%+.3f,%+.3f]\n\n', ...
            '', ci_b(1),ci_b(2), ternary(prod(ci_b)>0,'',' ⚠️含0'), ci_f(1),ci_f(2));
end


%% ================================================================
%  [C] 偏相关：控制 ROI 电极插值数（修正版 4c）
%% ================================================================
fprintf('\n===== [C] 偏相关（协变量 = 该被试 ROI 电极被插值数）=====\n');
QR = readtable('QC_RAT1_vs_RAT2.csv');
QT = readtable('QC_REST_vs_TRAIN.csv');
[okR,locR] = ismember(subs, QR.Subject);
[okT,locT] = ismember(subs, QT.Subject);
assert(all(okR) && all(okT), 'QC csv 缺被试');
covR = QR.nROIinterp(locR);      % RAT run
covT = QT.nROIinterp(locT);      % REST+TRAIN run
COV  = {covT; covT; covR; max(covT,covR)};   % 各检验对应的数据来源

for t = 1:size(TESTS,1)
    fprintf('  %-14s', TESTS{t,1});
    for g = 1:2
        k = (isFA == (g==2));
        x = TESTS{t,2}(k);  y = TESTS{t,3}(k);  c = COV{t}(k);
        [r0,p0] = corr(x,y);
        if std(c)==0
            fprintf('  %s: r=%+.3f (协变量无变异)', ternary(g==1,'CT','FA'), r0);
        else
            [r1,p1] = partialcorr(x,y,c);
            fprintf('  %s: r=%+.3f p=%.3f → 偏 r=%+.3f p=%.3f', ...
                    ternary(g==1,'CT','FA'), r0,p0, r1,p1);
        end
    end
    fprintf('\n');
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end