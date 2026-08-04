


%% Parse_Preproc_Summary.m
%  从两个 preprocessing_summary_*.txt 汇总预处理质控指标
%  输出：通道剔除统计、组间比较、先验 ROI 电极插值率、时长/epoch 数、排除名单
%  依赖：Statistics and Machine Learning Toolbox (ttest2, fishertest)

clear; clc;

% ===================== 配置 =====================
% {运行名, 文件路径, 状态1名, 状态2名, 状态1标称时长(s), 状态2标称时长(s)}
runs = { 'RAT1_vs_RAT2',  'preprocessing_summary_RAT1_RAT2.txt',  'RAT1', 'RAT2',  360, 360;
         'REST_vs_TRAIN', 'preprocessing_summary_REST_TRAIN.txt', 'REST', 'TRAIN', 120, 900 };

EPOCH_LEN = 20;                                   % wPLI 分段长度 (s)

% 先验假设涉及的电极
roi_H1 = {'Pz','F3','F4','Fp1','Fp2'};            % H1: Pz-前额 DMN
roi_H2 = {'F7','FC1','FC5'};                      % H2: 左额中语义-执行
roi_all = [roi_H1, roi_H2];

% 敏感性分析：想据此生成排除名单的电极（可改）
excl_targets = { 'RAT1_vs_RAT2',  {'F7'};         % H2 排除 F7 插值者
                 'REST_vs_TRAIN', {'Fp1','Fp2'} };

OUT_CSV = true;                                   % 是否导出逐被试 CSV

% ===================== 主循环 =====================
ALL = struct();

for r = 1:size(runs,1)
    runName = runs{r,1};  fname = runs{r,2};
    s1 = runs{r,3};       s2 = runs{r,4};
    d1_nom = runs{r,5};   d2_nom = runs{r,6};

    if ~exist(fname,'file')
        warning('未找到文件 %s，跳过。', fname); continue;
    end

    T = parse_summary(fname);
    n = height(T);
    isFA = startsWith(T.Subject,'FA');
    isCT = ~isFA;
    ALL.(matlab.lang.makeValidName(runName)) = T;

    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('  %s   (n = %d:  CT = %d,  FA = %d)\n', runName, n, sum(isCT), sum(isFA));
    fprintf('==========================================================\n');

    % ---------- 1. 通道剔除数 ----------
    k = T.nChanRemoved;
    fprintf('\n[1] 剔除/插值通道数\n');
    fprintf('    全体: M = %.2f, SD = %.2f, range %d–%d\n', mean(k), std(k), min(k), max(k));
    fprintf('    CT  : M = %.2f, SD = %.2f, range %d–%d\n', ...
            mean(k(isCT)), std(k(isCT)), min(k(isCT)), max(k(isCT)));
    fprintf('    FA  : M = %.2f, SD = %.2f, range %d–%d\n', ...
            mean(k(isFA)), std(k(isFA)), min(k(isFA)), max(k(isFA)));
    [~,p,~,st] = ttest2(k(isFA), k(isCT), 'Vartype','unequal');
    fprintf('    Welch t(%.1f) = %.2f, p = %.3f%s\n', st.df, st.tstat, p, star(p));
    fprintf('    >20%% 通道被剔除的被试 (31导 → >6.2): %s\n', ...
            listOrNone(T.Subject(k > 0.20*31)));

    % ---------- 2. ROI 电极插值率 ----------
    fprintf('\n[2] 先验 ROI 电极插值情况\n');
    fprintf('    %-6s %-4s %8s %8s %8s   %-9s  %s\n', ...
            '电极','假设','全体','CT','FA','Fisher p','被插值的被试');
    fprintf('    %s\n', repmat('-',1,96));
    roiHit = false(n, numel(roi_all));
    for c = 1:numel(roi_all)
        ch  = roi_all{c};
        hit = cellfun(@(x) any(strcmpi(x, ch)), T.RemovedList);
        roiHit(:,c) = hit;
        tbl = [sum(hit&isFA), sum(~hit&isFA); sum(hit&isCT), sum(~hit&isCT)];
        try, [~, pf] = fishertest(tbl); catch, pf = NaN; end
        hyp = 'H1'; if any(strcmp(ch, roi_H2)), hyp = 'H2'; end
        fprintf('    %-6s %-4s %3d(%4.1f%%) %3d(%4.1f%%) %3d(%4.1f%%)   %-9s  %s\n', ...
            ch, hyp, sum(hit), 100*mean(hit), ...
            sum(hit&isCT), 100*sum(hit&isCT)/sum(isCT), ...
            sum(hit&isFA), 100*sum(hit&isFA)/sum(isFA), ...
            sprintf('%.3f%s', pf, star(pf)), abbrevList(T.Subject(hit)));
    end
    nROI = sum(roiHit,2);
    fprintf('    %s\n', repmat('-',1,96));
    fprintf('    合计: %d / %d 人次 (%.1f%%)\n', sum(roiHit(:)), numel(roiHit), 100*mean(roiHit(:)));
    fprintf('    至少一个 ROI 电极被插值: %d 人 (CT %d, FA %d)\n', ...
            sum(nROI>0), sum(nROI>0 & isCT), sum(nROI>0 & isFA));
    [~,pn,~,stn] = ttest2(nROI(isFA), nROI(isCT), 'Vartype','unequal');
    fprintf('    每人 ROI 插值数: CT M = %.2f, FA M = %.2f; Welch t(%.1f) = %.2f, p = %.3f%s\n', ...
            mean(nROI(isCT)), mean(nROI(isFA)), stn.df, stn.tstat, pn, star(pn));
    T.nROIinterp = nROI;    % 可作偏相关协变量
    ALL.(matlab.lang.makeValidName(runName)) = T;

    % ---------- 3. 时长与 epoch 数 ----------
    e1 = floor(T.Dur1/EPOCH_LEN);  e2 = floor(T.Dur2/EPOCH_LEN);
    fprintf('\n[3] 清洗后时长与 %d s epoch 数\n', EPOCH_LEN);
    fprintf('    %-6s: %.1f–%.1f s  →  %d–%d epochs (标称 %d s = %d epochs)\n', ...
            s1, min(T.Dur1), max(T.Dur1), min(e1), max(e1), d1_nom, floor(d1_nom/EPOCH_LEN));
    fprintf('    %-6s: %.1f–%.1f s  →  %d–%d epochs (标称 %d s = %d epochs)\n', ...
            s2, min(T.Dur2), max(T.Dur2), min(e2), max(e2), d2_nom, floor(d2_nom/EPOCH_LEN));
    odd = T.Dur1 ~= d1_nom | T.Dur2 ~= d2_nom;
    if any(odd)
        fprintf('    时长偏离标称值的被试:\n');
        for i = find(odd)'
            fprintf('      %-6s  %s = %.0f s (%d ep),  %s = %.0f s (%d ep)\n', ...
                    T.Subject{i}, s1, T.Dur1(i), e1(i), s2, T.Dur2(i), e2(i));
        end
    else
        fprintf('    全部被试时长均为标称值 ✅\n');
    end

    % ---------- 4. 生成排除名单 ----------
    idx = find(strcmp(excl_targets(:,1), runName), 1);
    if ~isempty(idx)
        tg  = excl_targets{idx,2};
        hit = cellfun(@(x) any(ismember(lower(x), lower(tg))), T.RemovedList);
        fprintf('\n[4] 敏感性分析排除名单 —— 任一 {%s} 被插值\n', strjoin(tg,', '));
        fprintf('    排除 %d 人 (CT %d, FA %d) → 保留 CT %d, FA %d\n', ...
                sum(hit), sum(hit&isCT), sum(hit&isFA), sum(~hit&isCT), sum(~hit&isFA));
        fprintf('    可直接复制到分析脚本:\n');
        fprintf('    excl_list = {%s};\n', ...
                strjoin(cellfun(@(s) sprintf('''%s''',s), T.Subject(hit), 'uni',0), ', '));
    end

    % ---------- 5. 导出 ----------
    if OUT_CSV
        Tout = T(:, {'Subject','nChanRemoved','nROIinterp','Dur1','Dur2'});
        Tout.Group     = repmat({'CT'}, n, 1);  Tout.Group(isFA) = {'FA'};
        Tout.Removed   = cellfun(@(x) strjoin(x,'|'), T.RemovedList, 'uni', 0);
        Tout.nEpoch1   = e1;  Tout.nEpoch2 = e2;
        Tout.Properties.VariableNames{'Dur1'} = [s1 '_dur_s'];
        Tout.Properties.VariableNames{'Dur2'} = [s2 '_dur_s'];
        csv = sprintf('QC_%s.csv', runName);
        writetable(Tout, csv);
        fprintf('\n    📄 已导出: %s\n', csv);
    end
end

% ===================== 跨运行一致性检查 =====================
fn = fieldnames(ALL);
if numel(fn) == 2
    A = ALL.(fn{1});  B = ALL.(fn{2});
    fprintf('\n==========================================================\n');
    fprintf('  跨运行一致性检查 (Figure 3 依赖此项)\n');
    fprintf('==========================================================\n');
    if isequal(A.Subject, B.Subject)
        fprintf('  ✅ 两次运行的被试列表与顺序完全一致 (n = %d)\n', height(A));
    else
        fprintf('  ⚠️ 被试列表不一致！\n');
        fprintf('     仅在 %s: %s\n', fn{1}, listOrNone(setdiff(A.Subject, B.Subject)));
        fprintf('     仅在 %s: %s\n', fn{2}, listOrNone(setdiff(B.Subject, A.Subject)));
        fprintf('     ⚠️ 跨数据集分析前必须按 final_subs 显式对齐！\n');
    end
    fprintf('  提示: .mat 文件的顺序由 valid_subs_mask 决定，仍需在分析脚本中断言:\n');
    fprintf('        assert(isequal(A.final_subs, B.final_subs))\n');
end

fprintf('\n注: ICA 剔除成分数不在本 txt 中，需从 [subID]_preproc_report.mat 提取。\n\n');


%% ===================== 局部函数 =====================
function T = parse_summary(fname)
% 解析形如: "CT001: 移除通道[F7], EEG1=360.0s, EEG2=360.0s"
% 不匹配任何中文字符，因此不受 UTF-8 / GBK 编码差异影响
    fid = fopen(fname, 'r');
    if fid < 0, error('无法打开 %s', fname); end
    C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = C{1};

    sub = {}; rmv = {}; d1 = []; d2 = [];
    for i = 1:numel(lines)
        L = lines{i};
        if isempty(strtrim(L)), continue; end

        s = regexp(L, '([A-Za-z]{2}\d+)\s*:', 'tokens', 'once');
        if isempty(s), continue; end

        % 方括号内容：只保留 ASCII 电极名，"无" 自动被过滤为空
        b = regexp(L, '\[(.*?)\]', 'tokens', 'once');
        ch = {};
        if ~isempty(b)
            parts = strtrim(strsplit(b{1}, {',', '、', ';'}));
            ok = ~cellfun(@isempty, regexp(parts, '^[A-Za-z]+\d*$', 'once'));
            ch = unique(parts(ok), 'stable');
        end

        v1 = regexp(L, 'EEG1\s*=\s*([\d.]+)\s*s', 'tokens', 'once');
        v2 = regexp(L, 'EEG2\s*=\s*([\d.]+)\s*s', 'tokens', 'once');

        sub{end+1,1} = s{1};                                   %#ok<AGROW>
        rmv{end+1,1} = ch;                                     %#ok<AGROW>
        d1(end+1,1)  = ternary(isempty(v1), NaN, str2double(v1)); %#ok<AGROW>
        d2(end+1,1)  = ternary(isempty(v2), NaN, str2double(v2)); %#ok<AGROW>
    end

    if isempty(sub), error('未从 %s 解析到任何记录，请检查文件格式。', fname); end
    T = table(sub, rmv, cellfun(@numel, rmv), d1, d2, ...
        'VariableNames', {'Subject','RemovedList','nChanRemoved','Dur1','Dur2'});
end

function s = star(p)
    if isnan(p),      s = '';
    elseif p < .001,  s = ' ***';
    elseif p < .01,   s = ' **';
    elseif p < .05,   s = ' *';
    elseif p < .10,   s = ' †';
    else,             s = '';
    end
end

function s = listOrNone(c)
    if isempty(c), s = '(无)'; else, s = strjoin(c(:)', ', '); end
end

function s = abbrevList(c)
    if isempty(c), s = '—';
    elseif numel(c) <= 5, s = strjoin(c(:)', ',');
    else, s = sprintf('%s … (共%d人)', strjoin(c(1:4)', ','), numel(c));
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end

