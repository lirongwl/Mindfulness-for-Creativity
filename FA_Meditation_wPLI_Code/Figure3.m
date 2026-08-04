%% Plot_Scatter_Correlations_RealData_CrossState_Only
% 基于真实数据的 1x2 散点图 (仅展示跨状态：TRAIN-REST × RAT2-RAT1)
% 契合 JoCN 风格：精简冗余面板，去除图面交互作用 P 值，突出核心逻辑
%
% ★ 本版改动：
%   1. 仅保留原图的 b 和 e（即 FA 和 CT 的跨状态脑-脑相关）。
%   2. 布局改为 1x2，添加组别 Title。
%   3. 移除了图面上红色的 Interaction 标注（移至正文/图注中报告）。
%   4. 依然保留完整的统计计算和 Excel 导出，方便撰写正文。

clear; clc; close all;

% =========================================================================
% 0. 口径开关
% =========================================================================
OUTLIER_MODE = 'z3';        % 'raw' 主分析 | 'z3' 补充材料

switch OUTLIER_MODE
    case 'raw', global_z_thresh = Inf;   % abs(z)<=Inf 恒真 → 不剔除
    case 'z3',  global_z_thresh = 3;
    otherwise,  error('OUTLIER_MODE 只能是 ''raw'' 或 ''z3''');
end
fprintf('========================================================\n');
fprintf('  离群点口径: %s   (z 阈值 = %g)\n', OUTLIER_MODE, global_z_thresh);
fprintf('========================================================\n');

% =========================================================================
% 1. 行为学数据
% =========================================================================
CT_subs = {'CT001','CT004','CT016','CT022','CT034','CT035','CT042','CT055','CT061','CT071', ...
           'CT077','CT080','CT090','CT094','CT096','CT103','CT105','CT106','CT113','CT116'};
FA_subs = {'FA002','FA003','FA006','FA009','FA010','FA012','FA013','FA024','FA025','FA033', ...
           'FA050','FA059','FA062','FA063','FA070','FA075','FA078','FA092','FA098','FA101', ...
           'FA112','FA118','FA119'};

CT_RAT1 = [26, 23, 23, 17, 25, 26, 18, 23, 29, 31, 22, 24, 26, 19, 29, 29, 28, 27, 29, 25]';
CT_RAT2 = [33, 24, 25, 26, 25, 25, 30, 24, 23, 33, 23, 19, 31, 23, 24, 31, 29, 28, 33, 32]';
FA_RAT1 = [20, 27, 29, 24, 18, 31, 29, 13, 13, 23,  8, 24, 20, 21, 13, 27, 20, 25, 23, 23, 28, 22, 25]';
FA_RAT2 = [22, 32, 32, 30, 20, 33, 24, 24, 15, 35, 19, 23, 24, 20, 26, 32, 21, 25, 22, 23, 28, 25, 27]';

CT_DIFF_score = CT_RAT2 - CT_RAT1;
FA_DIFF_score = FA_RAT2 - FA_RAT1;

BEH = containers.Map([CT_subs, FA_subs], num2cell([CT_DIFF_score; FA_DIFF_score]));

% =========================================================================
% 2. 读取脑电数据 + 完整性断言
% =========================================================================
path_A = 'Individual_wPLI_Data_TrainRest_Alpha.mat';   % TRAIN - REST
path_B = 'Individual_wPLI_Data_RAT1RAT2_Alpha.mat';    % RAT2  - RAT1

if ~exist(path_A, 'file') || ~exist(path_B, 'file')
    error('未找到数据文件，请确保 .mat 文件在当前路径下！');
end

data_A = load(path_A);
data_B = load(path_B);

% ★ 关键断言：两个数据集必须同被试、同顺序、同通道
assert(isequal(data_A.final_subs, data_B.final_subs), ...
    ['❌ 两个数据集的被试顺序不一致！\n' ...
     '   必须先按 final_subs 显式对齐。']);
assert(isequal(data_A.chan_labels, data_B.chan_labels), '❌ 两个数据集通道顺序不一致！');

subs      = data_A.final_subs(:);
num_subs  = numel(subs);
chan_names = data_A.chan_labels;

% ★ 分组按被试名判定
isFA   = startsWith(subs, 'FA');
idx_CT = find(~isFA);
idx_FA = find(isFA);

% ★ 行为学按名查表
dRAT = nan(num_subs, 1);
for i = 1:num_subs
    assert(isKey(BEH, subs{i}), '❌ 被试 %s 不在行为学表中', subs{i});
    dRAT(i) = BEH(subs{i});
end

fprintf('✅ 数据完整性检查通过: n = %d (CT %d, FA %d)，脑–行为按名对齐\n\n', ...
        num_subs, numel(idx_CT), numel(idx_FA));

% 提取目标边
iF3  = find(strcmp(chan_names, 'F3'));
iPz  = find(strcmp(chan_names, 'Pz'));
iF7  = find(strcmp(chan_names, 'F7'));
iFC1 = find(strcmp(chan_names, 'FC1'));

edge_F3_Pz_all  = squeeze(data_A.final_wPLI_DIFF(:, iF3, iPz));    % TRAIN - REST
edge_F7_FC1_all = squeeze(data_B.final_wPLI_DIFF(:, iF7, iFC1));   % RAT2  - RAT1

CT_F3_Pz  = edge_F3_Pz_all(idx_CT);   FA_F3_Pz  = edge_F3_Pz_all(idx_FA);
CT_F7_FC1 = edge_F7_FC1_all(idx_CT);  FA_F7_FC1 = edge_F7_FC1_all(idx_FA);
CT_DIFF   = dRAT(idx_CT);             FA_DIFF   = dRAT(idx_FA);

% =========================================================================
% 3. 统计计算 (保留全部3组对比，用于控制台输出和Excel)
% =========================================================================
% Col 1: F3-Pz  vs  ΔRAT
[r_FA_1, p_FA_1, n_FA_1, ~, ~] = clean_corr_data(FA_F3_Pz, FA_DIFF, global_z_thresh);
[r_CT_1, p_CT_1, n_CT_1, ~, ~] = clean_corr_data(CT_F3_Pz, CT_DIFF, global_z_thresh);
[p_diff_1, z_diff_1, ci_1] = fisher_z_diff(r_CT_1, r_FA_1, n_CT_1, n_FA_1);

% Col 2: F3-Pz  vs  F7-FC1 (核心画图数据)
[r_FA_2, p_FA_2, n_FA_2, FA_F3_c2, FA_F7_c2] = clean_corr_data(FA_F3_Pz, FA_F7_FC1, global_z_thresh);
[r_CT_2, p_CT_2, n_CT_2, CT_F3_c2, CT_F7_c2] = clean_corr_data(CT_F3_Pz, CT_F7_FC1, global_z_thresh);
[p_diff_2, z_diff_2, ci_2] = fisher_z_diff(r_CT_2, r_FA_2, n_CT_2, n_FA_2);

% Col 3: F7-FC1 vs  ΔRAT
[r_FA_3, p_FA_3, n_FA_3, ~, ~] = clean_corr_data(FA_F7_FC1, FA_DIFF, global_z_thresh);
[r_CT_3, p_CT_3, n_CT_3, ~, ~] = clean_corr_data(CT_F7_FC1, CT_DIFF, global_z_thresh);
[p_diff_3, z_diff_3, ci_3] = fisher_z_diff(r_CT_3, r_FA_3, n_CT_3, n_FA_3);

% -------------------------------------------------------------------------
% 控制台汇总（可直接抄进稿子）
% -------------------------------------------------------------------------
col_names = {'F3-Pz vs dRAT', 'F3-Pz vs F7-FC1', 'F7-FC1 vs dRAT'};
RF = [r_FA_1 r_FA_2 r_FA_3];  PF = [p_FA_1 p_FA_2 p_FA_3];  NF = [n_FA_1 n_FA_2 n_FA_3];
RC = [r_CT_1 r_CT_2 r_CT_3];  PC = [p_CT_1 p_CT_2 p_CT_3];  NC = [n_CT_1 n_CT_2 n_CT_3];
PD = [p_diff_1 p_diff_2 p_diff_3];  ZD = [z_diff_1 z_diff_2 z_diff_3];
CI = [ci_1; ci_2; ci_3];

fprintf('---------------- Figure 3 统计汇总 (%s) ----------------\n', OUTLIER_MODE);
fprintf('%-18s | %-26s | %-26s | %s\n', '对比', 'FA', 'CT', 'Interaction');
fprintf('%s\n', repmat('-', 1, 104));
for c = 1:3
    fprintf('%-18s | r=%+.3f p=%.3f (n=%2d) | r=%+.3f p=%.3f (n=%2d) | z=%+.2f p=%.3f, dz 95%%CI [%+.2f,%+.2f]\n', ...
        col_names{c}, RF(c), PF(c), NF(c), RC(c), PC(c), NC(c), ZD(c), PD(c), CI(c,1), CI(c,2));
end
fprintf('%s\n\n', repmat('-', 1, 104));

% 导出统计表
StatsTable = table(col_names', RF', PF', NF', RC', PC', NC', ZD', PD', CI(:,1), CI(:,2), ...
    'VariableNames', {'Comparison','FA_r','FA_p','FA_n','CT_r','CT_p','CT_n', ...
                      'Fisher_z','Interaction_p','dz_CI_low','dz_CI_high'});
xls_name = sprintf('Between_Group_Statistics_Pz&FC1_Networks_%s.xlsx', OUTLIER_MODE);
writetable(StatsTable, xls_name);
fprintf('📊 统计表已导出: %s\n', xls_name);

% =========================================================================
% 4. 绘图排版 (精简为 1x2 布局)
% =========================================================================
color_scatter = [0.35, 0.55, 0.75];

% 调整 Figure 尺寸以适应 1x2 布局
fig = figure('Name', sprintf('Fig3_CrossState (%s)', OUTLIER_MODE), ...
             'Position', [50, 50, 1100, 500], 'Color', 'w');

w = 0.35; h = 0.70;
pos_A = [0.10, 0.18, w, h];   % 左图：FA
pos_B = [0.60, 0.18, w, h];   % 右图：CT

xlim_F3    = [min(edge_F3_Pz_all)-0.05,  max(edge_F3_Pz_all)+0.05];
xlim_F7    = [min(edge_F7_FC1_all)-0.05, max(edge_F7_FC1_all)+0.05];

label_F3    = '\Delta wPLI (F3 - Pz) [TRAIN-REST]';
label_F7    = '\Delta wPLI (F7 - FC1) [RAT2-RAT1]';

% =========================================================================
% 5. 绘制 Panels
% =========================================================================
% --- Panel a: FA Group ---
draw_scatter(pos_A, FA_F3_c2, FA_F7_c2, r_FA_2, p_FA_2, n_FA_2, label_F3, label_F7, color_scatter, xlim_F3, xlim_F7);
title('FA Group', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');

% --- Panel b: CT Group ---
draw_scatter(pos_B, CT_F3_c2, CT_F7_c2, r_CT_2, p_CT_2, n_CT_2, label_F3, label_F7, color_scatter, xlim_F3, xlim_F7);
title('CT Group', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');

% =========================================================================
% 6. 全局标签
% =========================================================================
axes('Position', [0 0 1 1], 'Visible', 'off');

label_props = {'FontSize', 28, 'FontWeight', 'bold', 'FontName', 'Arial'};
text(0.02, 0.95, 'a', label_props{:});   
text(0.52, 0.95, 'b', label_props{:});

% 口径水印（投稿前可注释掉）
text(0.99, 0.02, sprintf('[%s]', OUTLIER_MODE), 'FontSize', 8, 'Color', [0.6 0.6 0.6], ...
     'HorizontalAlignment', 'right');

% =========================================================================
% 7. 导出
% =========================================================================
pdf_name = sprintf('Fig3_CrossState_Only_%s.pdf', OUTLIER_MODE);
png_name = sprintf('Fig3_CrossState_Only_%s.png', OUTLIER_MODE);
exportgraphics(fig, pdf_name, 'ContentType', 'vector');
exportgraphics(fig, png_name, 'Resolution', 600);
fprintf('🎉 图片已保存: %s / %s\n', pdf_name, png_name);
fprintf('\n提示: 把 OUTLIER_MODE 改成 ''z3'' 再跑一次，即可得到补充材料版本。\n\n');


% =========================================================================
% 局部函数
% =========================================================================

function [r, p, n_valid, X_clean, Y_clean] = clean_corr_data(X, Y, z_thresh)
% z_thresh = Inf 时不剔除任何点（原始口径）
    X = X(:); Y = Y(:);
    if isinf(z_thresh)
        valid_idx = true(size(X));
    else
        Z_X = (X - mean(X)) / std(X);
        Z_Y = (Y - mean(Y)) / std(Y);
        valid_idx = abs(Z_X) <= z_thresh & abs(Z_Y) <= z_thresh;
    end
    X_clean = X(valid_idx);
    Y_clean = Y(valid_idx);
    n_valid = numel(X_clean);

    if n_valid > 2
        [r, p] = corr(X_clean, Y_clean, 'Type', 'Pearson');
    else
        r = NaN; p = NaN;
    end
end


function [p_diff, z_diff, ci] = fisher_z_diff(r1, r2, n1, n2)
% r1 = CT, r2 = FA；返回 p 值、z 统计量、Δz 的 95% CI (FA − CT)
    if isnan(r1) || isnan(r2) || n1 <= 3 || n2 <= 3
        p_diff = NaN; z_diff = NaN; ci = [NaN NaN]; return;
    end
    r1 = max(min(r1, 0.9999), -0.9999);
    r2 = max(min(r2, 0.9999), -0.9999);
    z1 = atanh(r1);   % CT
    z2 = atanh(r2);   % FA
    SE_diff = sqrt(1/(n1-3) + 1/(n2-3));
    d       = z2 - z1;                  % FA − CT
    z_diff  = d / SE_diff;
    p_diff  = 2 * (1 - normcdf(abs(z_diff)));
    ci      = d + [-1 1] * 1.96 * SE_diff;
end


function draw_scatter(ax_pos, x_data, y_data, r_val, p_val, n_val, ...
                      x_label_str, y_label_str, color_scatter, x_lims, y_lims)
    axes('Position', ax_pos); hold on;

    scatter(x_data, y_data, 65, color_scatter, 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w', 'LineWidth', 1);

    if numel(x_data) > 2 && ~isnan(r_val)
        p_fit = polyfit(x_data, y_data, 1);
        x_fit = linspace(x_lims(1), x_lims(2), 100);
        y_fit = polyval(p_fit, x_fit);
        y_res = y_data - polyval(p_fit, x_data);
        n     = numel(x_data);
        SE    = sqrt(sum(y_res.^2) / (n - 2));
        tcrit = tinv(0.975, n - 2);          % 小样本用 t 分位数，不用 1.96
        CI_off = tcrit * SE * sqrt(1/n + (x_fit - mean(x_data)).^2 / sum((x_data - mean(x_data)).^2));

        patch([x_fit, fliplr(x_fit)], [y_fit + CI_off, fliplr(y_fit - CI_off)], ...
              [0.6 0.6 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(x_fit, y_fit, 'Color', [0.15 0.15 0.15], 'LineWidth', 2.5);
    end

    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, ...
             'FontName', 'Arial', 'FontSize', 13);
    xlim(x_lims); ylim(y_lims);
    xlabel(x_label_str, 'FontSize', 15, 'FontWeight', 'bold');
    ylabel(y_label_str, 'FontSize', 15, 'FontWeight', 'bold');

    if     p_val < 0.001, p_str = '***{\it P} < 0.001';
    elseif p_val < 0.01,  p_str = sprintf('**{\\it P} = %.3f', p_val);
    elseif p_val < 0.05,  p_str = sprintf('*{\\it P} = %.3f',  p_val);
    else,                 p_str = sprintf('{\\it P} = %.3f',   p_val);
    end

    % ★ 显示 n
    stats_text = sprintf('{\\it r} = %.3f\n%s', r_val, p_str);

    text(0.95, 0.92, stats_text, 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1, 'Margin', 6);

    hold off;
end