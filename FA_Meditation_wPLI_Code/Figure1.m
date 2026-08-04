%% Plot_Network_Correlations_ROI_.m

clear; clc; close all;

% =========================================================================
% 0. 口径开关
% =========================================================================
OUTLIER_MODE = 'raw';        % 'raw' 主分析 | 'z3' 补充材料

switch OUTLIER_MODE
    case 'raw', global_z_thresh = Inf;   % abs(z)<=Inf 恒真 → 不剔除
    case 'z3',  global_z_thresh = 3;
    otherwise,  error('OUTLIER_MODE 只能是 ''raw'' 或 ''z3''');
end
fprintf('========================================================\n');
fprintf('  离群点口径: %s   (z 阈值 = %g)\n', OUTLIER_MODE, global_z_thresh);
fprintf('========================================================\n');

% =========================================================================
% 1. 参数设置与数据准备 (使用真实数据)
% =========================================================================
dataFile = 'Individual_wPLI_Data_TrainRest_Alpha.mat';
if ~exist(dataFile, 'file')
    error('未找到数据文件 %s，请检查路径！', dataFile);
end
load(dataFile, 'final_wPLI_DIFF', 'chan_labels', 'idx_CT', 'idx_FA', 'final_subs');

% --- 录入真实的双盲行为学分数 (RAT1 和 RAT2) ---
CT_RAT1 = [26, 23, 23, 17, 25, 26, 18, 23, 29, 31, 22, 24, 26, 19, 29, 29, 28, 27, 29, 25]';
CT_RAT2 = [33, 24, 25, 26, 25, 25, 30, 24, 23, 33, 23, 19, 31, 23, 24, 31, 29, 28, 33, 32]';
CT_DIFF = CT_RAT2 - CT_RAT1;

FA_RAT1 = [20, 27, 29, 24, 18, 31, 29, 13, 13, 23, 8, 24, 20, 21, 13, 27, 20, 25, 23, 23, 28, 22, 25]';
FA_RAT2 = [22, 32, 32, 30, 20, 33, 24, 24, 15, 35, 19, 23, 24, 20, 26, 32, 21, 25, 22, 23, 28, 25, 27]';
FA_DIFF = FA_RAT2 - FA_RAT1;

% 目标电极
seed_node = 'Pz';
target_nodes = {'Fp1', 'Fp2', 'F3', 'F4'};
num_targets = length(target_nodes);

% 获取电极索引
idx_seed = find(strcmp(chan_labels, seed_node));
idx_targs = zeros(1, num_targets);
for i = 1:num_targets
    idx_targs(i) = find(strcmp(chan_labels, target_nodes{i}));
end

% 提取核心连接变化数据 (ΔwPLI)
wPLI_DIFF_CT = final_wPLI_DIFF(idx_CT, idx_seed, idx_targs);
wPLI_DIFF_FA = final_wPLI_DIFF(idx_FA, idx_seed, idx_targs);

% 计算全脑平均连接变化矩阵 (热图仅作描述性展示，保留全样本)
mean_matrix_CT = squeeze(mean(final_wPLI_DIFF(idx_CT, :, :), 1));
mean_matrix_FA = squeeze(mean(final_wPLI_DIFF(idx_FA, :, :), 1));

% =========================================================================
% 2. 统计计算 (包含组内和组间差异，已加入异常值剔除)
% =========================================================================
stats_CT = struct('mean_diff', zeros(1,4), 'r_corr', zeros(1,4), 'p_corr', zeros(1,4));
stats_FA = struct('mean_diff', zeros(1,4), 'r_corr', zeros(1,4), 'p_corr', zeros(1,4));

% 用于存储组间统计结果的数组
p_diff_wPLI_edges = zeros(1,4);
p_diff_corr_edges = zeros(1,4);
z_diff_corr_edges = zeros(1,4);

for i = 1:num_targets
    % CT 组计算 (根据口径剔除异常值)
    data_CT = squeeze(wPLI_DIFF_CT(:, 1, i));
    [data_CT_clean, CT_DIFF_clean, ~] = remove_outliers(data_CT, CT_DIFF, global_z_thresh);
    
    stats_CT.mean_diff(i) = mean(data_CT_clean);
    [r_CT_mat, p_CT_mat] = corrcoef(data_CT_clean, CT_DIFF_clean);
    stats_CT.r_corr(i) = r_CT_mat(1,2);
    stats_CT.p_corr(i) = p_CT_mat(1,2);
    
    % FA 组计算 (根据口径剔除异常值)
    data_FA = squeeze(wPLI_DIFF_FA(:, 1, i));
    [data_FA_clean, FA_DIFF_clean, ~] = remove_outliers(data_FA, FA_DIFF, global_z_thresh);
    
    stats_FA.mean_diff(i) = mean(data_FA_clean);
    [r_FA_mat, p_FA_mat] = corrcoef(data_FA_clean, FA_DIFF_clean);
    stats_FA.r_corr(i) = r_FA_mat(1,2);
    stats_FA.p_corr(i) = p_FA_mat(1,2);
    
    % --- 组间差异：连接变化量 (t-test) ---
    [~, p_diff_wPLI_edges(i)] = ttest2(data_FA_clean, data_CT_clean);
    
    % --- 组间差异：相关系数 (Fisher's Z transform) ---
    z_FA = 0.5 * log((1 + stats_FA.r_corr(i)) / (1 - stats_FA.r_corr(i)));
    z_CT = 0.5 * log((1 + stats_CT.r_corr(i)) / (1 - stats_CT.r_corr(i)));
    
    % 动态获取剔除异常值后的样本量
    n_FA_clean = length(data_FA_clean);
    n_CT_clean = length(data_CT_clean);
    
    SE_diff = sqrt((1 / (n_FA_clean - 3)) + (1 / (n_CT_clean - 3)));
    z_diff_corr_edges(i) = (z_FA - z_CT) / SE_diff;
    p_diff_corr_edges(i) = 2 * (1 - normcdf(abs(z_diff_corr_edges(i)))); % 双尾 p 值
end

% ROI 级别计算
wPLI_DIFF_CT_ROI = mean(squeeze(wPLI_DIFF_CT), 2);
wPLI_DIFF_FA_ROI = mean(squeeze(wPLI_DIFF_FA), 2);

% ROI 级别异常值剔除与相关性计算
[wPLI_CT_ROI_clean, CT_DIFF_ROI_clean, ~] = remove_outliers(wPLI_DIFF_CT_ROI, CT_DIFF, global_z_thresh);
[r_CT_ROI_mat, p_CT_ROI_mat] = corrcoef(wPLI_CT_ROI_clean, CT_DIFF_ROI_clean);
r_CT_ROI = r_CT_ROI_mat(1,2); p_CT_ROI = p_CT_ROI_mat(1,2);

[wPLI_FA_ROI_clean, FA_DIFF_ROI_clean, ~] = remove_outliers(wPLI_DIFF_FA_ROI, FA_DIFF, global_z_thresh);
[r_FA_ROI_mat, p_FA_ROI_mat] = corrcoef(wPLI_FA_ROI_clean, FA_DIFF_ROI_clean);
r_FA_ROI = r_FA_ROI_mat(1,2); p_FA_ROI = p_FA_ROI_mat(1,2);

% --- ROI 级别组间差异 ---
[~, p_diff_wPLI_ROI] = ttest2(wPLI_FA_ROI_clean, wPLI_CT_ROI_clean);

z_FA_ROI = 0.5 * log((1 + r_FA_ROI) / (1 - r_FA_ROI));
z_CT_ROI = 0.5 * log((1 + r_CT_ROI) / (1 - r_CT_ROI));
n_FA_ROI_clean = length(wPLI_FA_ROI_clean);
n_CT_ROI_clean = length(wPLI_CT_ROI_clean);
z_diff_corr_ROI = (z_FA_ROI - z_CT_ROI) / sqrt((1 / (n_FA_ROI_clean - 3)) + (1 / (n_CT_ROI_clean - 3)));
p_diff_corr_ROI = 2 * (1 - normcdf(abs(z_diff_corr_ROI)));

% =========================================================================
% 3. 导出统计结果到 Excel
% =========================================================================
Connection_Names = {'Pz-Fp1'; 'Pz-Fp2'; 'Pz-F3'; 'Pz-F4'; 'Pz-Prefrontal_ROI'};

FA_Mean_wPLI = [stats_FA.mean_diff, mean(wPLI_FA_ROI_clean)]';
CT_Mean_wPLI = [stats_CT.mean_diff, mean(wPLI_CT_ROI_clean)]';
P_Value_wPLI_Diff = [p_diff_wPLI_edges, p_diff_wPLI_ROI]';

FA_Correlation_r = [stats_FA.r_corr, r_FA_ROI]';
FA_Correlation_p = [stats_FA.p_corr, p_FA_ROI]';
CT_Correlation_r = [stats_CT.r_corr, r_CT_ROI]';
CT_Correlation_p = [stats_CT.p_corr, p_CT_ROI]';
Z_Value_Corr_Diff = [z_diff_corr_edges, z_diff_corr_ROI]';
P_Value_Corr_Diff = [p_diff_corr_edges, p_diff_corr_ROI]';

% 创建表格
StatsTable = table(Connection_Names, ...
    FA_Mean_wPLI, CT_Mean_wPLI, P_Value_wPLI_Diff, ...
    FA_Correlation_r, FA_Correlation_p, CT_Correlation_r, CT_Correlation_p, ...
    Z_Value_Corr_Diff, P_Value_Corr_Diff);

% 写入 Excel (文件名带上口径后缀)
excel_filename = sprintf('Between_Group_Statistics_Pz_Network_%s.xlsx', OUTLIER_MODE);
writetable(StatsTable, excel_filename);
fprintf('📊 组间统计结果已成功导出至: %s\n', excel_filename);

% 全局统一线宽和热图颜色映射标准
global_max_diff = max(max(abs(stats_CT.mean_diff)), max(abs(stats_FA.mean_diff)));
if global_max_diff == 0, global_max_diff = 1; end 

all_vals = [mean_matrix_CT(:); mean_matrix_FA(:)];
heatmap_cmax = prctile(abs(all_vals), 98); 
if heatmap_cmax == 0, heatmap_cmax = 0.1; end
heatmap_clims = [-heatmap_cmax, heatmap_cmax];

% =========================================================================
% 4. 全局视觉与排版设置
% =========================================================================
color_inc = [0.85, 0.20, 0.20]; 
color_dec = [0.15, 0.15, 0.15]; 
color_nonsig = [0.15, 0.15, 0.15]; 
color_scatter = [0.35, 0.55, 0.75]; % 调整为更接近附图2的灰蓝色

% 加宽画布，为左侧的全局行标签留出空间
fig_name = sprintf('Network & Correlation (%s)', OUTLIER_MODE);
fig = figure('Name', fig_name, 'Position', [50, 50, 1650, 850], 'Color', 'w');

% 重新分配 3x2 子图位置 [left, bottom, width, height]
pos_A_heat  = [0.08, 0.55, 0.22, 0.35];
pos_B_brain = [0.38, 0.55, 0.23, 0.38];
pos_C_scat  = [0.70, 0.58, 0.25, 0.32];

pos_D_heat  = [0.08, 0.08, 0.22, 0.35];
pos_E_brain = [0.38, 0.08, 0.23, 0.38];
pos_F_scat  = [0.70, 0.11, 0.25, 0.32];

x_limits = [-0.6, 0.4];
y_limits = [-15, 20];

% =========================================================================
% 5. 组合绘制六大 Panels (仅第一行有标题)
% =========================================================================

% --- Row 1: FA Group ---
draw_heatmap(pos_A_heat, mean_matrix_FA, chan_labels, '', heatmap_clims);
draw_network(pos_B_brain, stats_FA, '', color_inc, color_dec, color_nonsig, global_max_diff);
% 注意：散点图传入的是剔除异常值后的 clean 数据
draw_scatter(pos_C_scat, wPLI_FA_ROI_clean, FA_DIFF_ROI_clean, r_FA_ROI, p_FA_ROI, 'Prefrontal ROI', '', color_scatter, x_limits, y_limits);

% --- Row 2: CT Group (标题留空) ---
draw_heatmap(pos_D_heat, mean_matrix_CT, chan_labels, '', heatmap_clims);
draw_network(pos_E_brain, stats_CT, '', color_inc, color_dec, color_nonsig, global_max_diff);
% 注意：散点图传入的是剔除异常值后的 clean 数据
draw_scatter(pos_F_scat, wPLI_CT_ROI_clean, CT_DIFF_ROI_clean, r_CT_ROI, p_CT_ROI, 'Prefrontal ROI', '', color_scatter, x_limits, y_limits);

% =========================================================================
% 6. 绘制全局行标签 (Row Labels) 与 a-f 字母标签
% =========================================================================
% 创建一个覆盖全图的透明坐标轴，用于放置全局文字
axes('Position', [0 0 1 1], 'Visible', 'off');

% 绘制左侧垂直的组别标签 (Row Labels)
text(0.02, 0.725, 'FA Group', 'Rotation', 90, 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial', 'HorizontalAlignment', 'center');
text(0.02, 0.255, 'CT Group', 'Rotation', 90, 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial', 'HorizontalAlignment', 'center');

% 绘制 a-f 标签
label_props = {'FontSize', 24, 'FontWeight', 'bold', 'FontName', 'Arial'};
text(0.03, 0.95, 'a', label_props{:});
text(0.34, 0.95, 'b', label_props{:});
text(0.65, 0.95, 'c', label_props{:});
text(0.03, 0.48, 'd', label_props{:});
text(0.34, 0.48, 'e', label_props{:});
text(0.65, 0.48, 'f', label_props{:});

% 在 c (FA组散点) 和 f (CT组散点) 之间添加组间差异显著性标注
text(0.825, 0.47, sprintf('Interaction: {\\it P} = %.3f', p_diff_corr_ROI), ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', 'Color', [0.8 0.1 0.1], ... % 红色字体
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', 'w');

% 口径水印（投稿前可注释掉）
text(0.99, 0.02, sprintf('[%s]', OUTLIER_MODE), 'FontSize', 8, 'Color', [0.6 0.6 0.6], ...
     'HorizontalAlignment', 'right');

% 保存高分辨率图片 (文件名带上口径后缀)
fig_save_name = sprintf('Network_Correlations_Pz_%s', OUTLIER_MODE);
savefig(fig, [fig_save_name, '.fig']);
exportgraphics(fig, [fig_save_name, '.png'], 'Resolution', 600);
fprintf('🎉 顶刊排版级高质量图片已保存: %s.png\n', fig_save_name);

% =========================================================================
% 局部函数定义区
% =========================================================================

function [x_clean, y_clean, keep_idx] = remove_outliers(x, y, z_thresh)
    % 基于 z_thresh 法则剔除双变量异常值
    % 如果 z_thresh 是 Inf，则保留全部数据
    if isinf(z_thresh)
        keep_idx = true(size(x));
    else
        zx = abs((x - mean(x)) / std(x));
        zy = abs((y - mean(y)) / std(y));
        keep_idx = (zx <= z_thresh) & (zy <= z_thresh);
    end
    x_clean = x(keep_idx);
    y_clean = y(keep_idx);
end

function draw_heatmap(ax_pos, data_matrix, labels, title_str, clims)
    axes('Position', ax_pos);
    imagesc(data_matrix);
    colormap(jet); 
    caxis(clims);  
    
    num_chans = length(labels);
    xticks(1:num_chans); yticks(1:num_chans);
    xticklabels(labels); yticklabels(labels);
    xtickangle(90);
    set(gca, 'TickDir', 'out', 'FontSize', 8, 'FontName', 'Arial');
    
    cb = colorbar;
    cb.LineWidth = 1; cb.FontSize = 10;
    
    if ~isempty(title_str)
        title(title_str, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
end

function draw_network(ax_pos, stats, title_str, color_inc, color_dec, color_nonsig, global_max_diff)
    axes('Position', ax_pos); hold on;
    nodes_info = {'Pz', 0.00, -0.55; 'Fp1', -0.22, 0.80; 'Fp2', 0.22, 0.80; 'F3', -0.42, 0.52; 'F4', 0.42, 0.52};
    
    if exist('brain.png', 'file')
        img = imread('brain.png');
        gray_img = rgb2gray(img); 
        mask = gray_img < 20; 
        r = img(:,:,1); r(mask) = 255;
        g = img(:,:,2); g(mask) = 255;
        b = img(:,:,3); b(mask) = 255;
        clean_img = cat(3, r, g, b);
        
        [rows, cols] = find(~mask);
        min_y = min(rows); max_y = max(rows); min_x = min(cols); max_x = max(cols);
        cropped_img = clean_img(min_y:max_y, min_x:max_x, :);
        imshow(cropped_img); hold on;
        
        brain_w = max_x - min_x + 1; brain_h = max_y - min_y + 1;
        brain_center_x = brain_w / 2; brain_center_y = brain_h / 2;
        node_x = brain_center_x + (cell2mat(nodes_info(:,2)) .* (brain_w / 2));
        node_y = brain_center_y - (cell2mat(nodes_info(:,3)) .* (brain_h / 2 * 1.05));
    else
        set(gca, 'XLim', [-1 1], 'YLim', [-1 1], 'Color', 'w'); axis equal off;
        node_x = cell2mat(nodes_info(:,2)); node_y = cell2mat(nodes_info(:,3));
    end
    
    for i = 1:4
        m_diff = stats.mean_diff(i); p_val = stats.p_corr(i);
        l_width = 1.0 + (5.0 * (abs(m_diff) / global_max_diff)); 
        
        if p_val < 0.05
            if m_diff > 0, e_color = color_inc; else, e_color = color_dec; end
            plot([node_x(1), node_x(i+1)], [node_y(1), node_y(i+1)], 'Color', e_color, 'LineStyle', '-', 'LineWidth', l_width);
        else
            e_color = color_nonsig; 
            draw_fine_dashed_line(node_x(1), node_x(i+1), node_y(1), node_y(i+1), e_color, l_width);
        end
    end
    
    final_node_size = 431; 
    scatter(node_x(1), node_y(1), final_node_size, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2); 
    scatter(node_x(2:end), node_y(2:end), final_node_size, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2); 
    
    for i = 1:5
        text(node_x(i), node_y(i), nodes_info{i,1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Arial', 'FontSize', 9.5, 'FontWeight', 'bold');
    end
    
    if ~isempty(title_str)
        title(title_str, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    hold off;
end

function draw_fine_dashed_line(x1, x2, y1, y2, color, l_width)
    dash_len = 150; gap_len = 150;  
    dist = sqrt((x2-x1)^2 + (y2-y1)^2);
    num_dashes = floor(dist / (dash_len + gap_len));
    
    if num_dashes == 0
        plot([x1, x2], [y1, y2], 'Color', color, 'LineWidth', l_width); return;
    end
    
    dx = (x2 - x1) / dist; dy = (y2 - y1) / dist;
    curr_x = x1; curr_y = y1;
    
    for k = 1:num_dashes
        end_x = curr_x + dash_len * dx; end_y = curr_y + dash_len * dy;
        plot([curr_x, end_x], [curr_y, end_y], 'Color', color, 'LineWidth', l_width, 'LineStyle', '-');
        curr_x = end_x + gap_len * dx; curr_y = end_y + gap_len * dy;
    end
    
    rem_dist = sqrt((x2-curr_x)^2 + (y2-curr_y)^2);
    if rem_dist > 0
        if rem_dist <= dash_len
            plot([curr_x, x2], [curr_y, y2], 'Color', color, 'LineWidth', l_width, 'LineStyle', '-');
        else
            plot([curr_x, curr_x + dash_len * dx], [curr_y, curr_y + dash_len * dy], 'Color', color, 'LineWidth', l_width, 'LineStyle', '-');
        end
    end
end

function draw_scatter(ax_pos, x_data, y_data, r_val, p_val, target_name, title_str, color_scatter, x_lims, y_lims)
    axes('Position', ax_pos); hold on;
    scatter(x_data, y_data, 65, color_scatter, 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'w', 'LineWidth', 1);
    
    p_fit = polyfit(x_data, y_data, 1);
    x_fit = linspace(x_lims(1), x_lims(2), 100);
    y_fit = polyval(p_fit, x_fit);
    y_resid = y_data - polyval(p_fit, x_data);
    SSresid = sum(y_resid.^2); 
    SE = sqrt(SSresid / (length(y_data)-2));
    CI_offset = 1.96 * SE * sqrt(1/length(x_data) + (x_fit - mean(x_data)).^2 / sum((x_data - mean(x_data)).^2));
    
    fill_color = [0.6 0.6 0.6]; line_color = [0.15 0.15 0.15]; 
    
    patch([x_fit, fliplr(x_fit)], [y_fit + CI_offset, fliplr(y_fit - CI_offset)], fill_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(x_fit, y_fit, 'Color', line_color, 'LineWidth', 2.5);
    
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 12);
    xlim(x_lims); ylim(y_lims);
    
    xlabel(sprintf('\\Delta wPLI (Pz - %s)', target_name), 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('\Delta RAT Score', 'FontSize', 15, 'FontWeight', 'bold');
    
    if ~isempty(title_str)
        title(title_str, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    
    if p_val < 0.001
        p_str = sprintf('***{\\it P} = %.4f', p_val);
    elseif p_val < 0.01
        p_str = sprintf('**{\\it P} = %.3f', p_val);
    elseif p_val < 0.05
        p_str = sprintf('*{\\it P} = %.3f', p_val);
    else
        p_str = sprintf('{\\it P} = %.3f', p_val);
    end   
    
    % 使用 TeX 语法将 r 和 P 变为斜体
    stats_text = sprintf('{\\it r} = %.3f\n%s', r_val, p_str);
    
    % 绘制带黑框白底的文本，绝对锚定在右上角 (x=0.95, y=0.92)
    text(0.95, 0.92, stats_text, 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1, 'Margin', 6);
        
    hold off;
end