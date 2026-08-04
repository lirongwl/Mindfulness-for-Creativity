%% Robustness_Checks.m
% 执行论文中描述的 4 种稳健性检验 (仅针对 FA 组)
% 1. Leave-one-out
% 2. Exclusion of participants with extensive channel reconstruction (>20%)
% 3. Rank-based estimation (Spearman)
% 4. Control for signal reconstruction (Partial correlation)

clear; clc; close all;

fprintf('========================================================\n');
fprintf('        Robustness Checks for FA Group Associations\n');
fprintf('========================================================\n\n');

% =========================================================================
% 1. 核心数据准备 (提取 FA 组)
% =========================================================================
% 读取两个状态的数据
data_TR = load('Individual_wPLI_Data_TrainRest_Alpha.mat'); % TRAIN-REST
data_RR = load('Individual_wPLI_Data_RAT1RAT2_Alpha.mat');  % RAT2-RAT1

% 真实的 FA 组被试列表 (共 23 人，与 TimeRecorder.m 和预处理记录一致)
FA_subs = {'FA002', 'FA003', 'FA006', 'FA009', 'FA010', 'FA012', 'FA013', ...
           'FA024', 'FA025', 'FA033', 'FA050', 'FA059', 'FA062', 'FA063', ...
           'FA070', 'FA075', 'FA078', 'FA092', 'FA098', 'FA101', 'FA112', ...
           'FA118', 'FA119'}';

% 行为学分数 (FA组, 23人)
FA_RAT1 = [20, 27, 29, 24, 18, 31, 29, 13, 13, 23, 8, 24, 20, 21, 13, 27, 20, 25, 23, 23, 28, 22, 25]';
FA_RAT2 = [22, 32, 32, 30, 20, 33, 24, 24, 15, 35, 19, 23, 24, 20, 26, 32, 21, 25, 22, 23, 28, 25, 27]';
dRAT_FA = FA_RAT2 - FA_RAT1;
n_FA = length(dRAT_FA);
idx_FA = data_TR.idx_FA;

% 提取 H1: Pz - Prefrontal ROI (TRAIN-REST)
idx_Pz = find(strcmp(data_TR.chan_labels, 'Pz'));
idx_Fp1 = find(strcmp(data_TR.chan_labels, 'Fp1'));
idx_Fp2 = find(strcmp(data_TR.chan_labels, 'Fp2'));
idx_F3 = find(strcmp(data_TR.chan_labels, 'F3'));
idx_F4 = find(strcmp(data_TR.chan_labels, 'F4'));
wPLI_Pz_ROI = mean(data_TR.final_wPLI_DIFF(idx_FA, idx_Pz, [idx_Fp1, idx_Fp2, idx_F3, idx_F4]), 3);
H1_ROI = squeeze(wPLI_Pz_ROI);

% 提取 H1: Pz - F3 (TRAIN-REST)
H1_PzF3 = squeeze(data_TR.final_wPLI_DIFF(idx_FA, idx_Pz, idx_F3));

% 提取 H2: FC1 - F7 (RAT2-RAT1)
idx_FC1 = find(strcmp(data_RR.chan_labels, 'FC1'));
idx_F7 = find(strcmp(data_RR.chan_labels, 'F7'));
H2_F7FC1 = squeeze(data_RR.final_wPLI_DIFF(idx_FA, idx_FC1, idx_F7));

% 提取 H3: Cross-state (F3-Pz vs F7-FC1)
H3_X = H1_PzF3;
H3_Y = H2_F7FC1;

% 汇总到一个结构体方便循环
vars = {'H1_ROI', 'H1_PzF3', 'H2_F7FC1', 'H3_Cross'};
X_data = {H1_ROI, H1_PzF3, H2_F7FC1, H3_X};
Y_data = {dRAT_FA, dRAT_FA, dRAT_FA, H3_Y};

% 计算基线 Pearson 相关系数 (用于对比)
base_r = zeros(1, 4);
for i = 1:4
    base_r(i) = corr(X_data{i}, Y_data{i}, 'Type', 'Pearson');
end

% =========================================================================
% 2. 导入真实的插值电极数据 (基于 preprocessing_summary 文本)
% =========================================================================
% 提取自 RAT1_RAT2 文件的移除通道数 (对应 FA_subs 顺序)
interp_RR = [2; 2; 7; 5; 5; 6; 7; 5; 2; 3; 5; 5; 5; 2; 7; 5; 11; 6; 4; 0; 2; 6; 2]; 

% 提取自 REST_TRAIN 文件的移除通道数 (对应 FA_subs 顺序)
interp_TR = [5; 1; 4; 2; 6; 3; 3; 1; 5; 2; 1; 2; 3; 1; 2; 1; 5; 0; 4; 1; 3; 4; 2];

% =========================================================================
% 3. 执行 Check 1: Leave-one-out
% =========================================================================
fprintf('--- Check 1: Leave-one-out ---\n');
for i = 1:4
    X = X_data{i}; Y = Y_data{i};
    r_loo = zeros(n_FA, 1);
    for skip = 1:n_FA
        keep = setdiff(1:n_FA, skip);
        r_loo(skip) = corr(X(keep), Y(keep), 'Type', 'Pearson');
    end
    min_r = min(r_loo); max_r = max(r_loo);
    max_delta = max(abs(r_loo - base_r(i)));
    fprintf('%-10s: r in [%+.2f, %+.2f], max Δr = %.2f\n', vars{i}, min_r, max_r, max_delta);
end
fprintf('\n');

% =========================================================================
% 4. 执行 Check 2: Exclusion of participants with extensive channel reconstruction
% =========================================================================
fprintf('--- Check 2: Exclude >20%% interpolated channels ---\n');
% 假设32导联，20%为6.4。因此剔除通道数 >= 7 的被试
% 我们以 RAT1_RAT2 阶段的插值数量作为剔除标准 (因为这4个人主要在这个阶段插值过多)
valid_idx = interp_RR <= 6; 
n_valid = sum(valid_idx);
excluded_subs = FA_subs(~valid_idx);

fprintf('Excluded %d participants: %s\n', length(excluded_subs), strjoin(excluded_subs, ', '));
fprintf('Retained %d/%d participants.\n', n_valid, n_FA);

for i = 1:4
    X = X_data{i}; Y = Y_data{i};
    X_clean = X(valid_idx); Y_clean = Y(valid_idx);
    [r, p] = corr(X_clean, Y_clean, 'Type', 'Pearson');
    fprintf('%-10s: r = %+.3f, p = %.3f (n = %d)\n', vars{i}, r, p, n_valid);
end
fprintf('\n');

% =========================================================================
% 5. 执行 Check 3: Rank-based estimation (Spearman)
% =========================================================================
fprintf('--- Check 3: Rank-based estimation (Spearman) ---\n');
for i = 1:4
    X = X_data{i}; Y = Y_data{i};
    [rho, p] = corr(X, Y, 'Type', 'Spearman');
    fprintf('%-10s: rho = %+.3f, p = %.3f\n', vars{i}, rho, p);
end
fprintf('\n');

% =========================================================================
% 6. 执行 Check 4: Control for signal reconstruction (Partial correlation)
% =========================================================================
fprintf('--- Check 4: Control for signal reconstruction ---\n');
for i = 1:4
    X = X_data{i}; Y = Y_data{i};
    
    % 根据假设所处的状态，选择对应的插值电极数量作为协变量
    if i == 1 || i == 2
        covariate = interp_TR; % H1 属于 TRAIN-REST 状态
    elseif i == 3
        covariate = interp_RR; % H2 属于 RAT2-RAT1 状态
    else
        covariate = max(interp_TR, interp_RR); % H3 跨状态，取两次中插值较多的一次作为最严苛控制
    end
    
    % 使用 partialcorr 控制插值电极数量
    [r_part, p_part] = partialcorr(X, Y, covariate);
    delta_r = r_part - base_r(i);
    fprintf('%-10s: partial r = %+.3f, p = %.3f (Δr = %+.3f)\n', vars{i}, r_part, p_part, delta_r);
end
fprintf('========================================================\n');


% =========================================================================
% 7. 执行 Check 5: Control for Sex (Partial correlation)
% =========================================================================
fprintf('\n--- Check 5: Control for Sex (Female=0, Male=1) ---\n');

% 录入 FA 组的性别数据 (严格按照 FA_subs 的顺序)
% FA_subs = {'FA002','FA003','FA006','FA009','FA010','FA012','FA013','FA024','FA025','FA033',...
%            'FA050','FA059','FA062','FA063','FA070','FA075','FA078','FA092','FA098','FA101',...
%            'FA112','FA118','FA119'};
% 男=1 (006, 013, 024, 063, 119), 女=0
FA_sex = [0, 0, 1, 0, 0, 0, 1, 1, 0, 0, ...
          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, ...
          0, 0, 1]';

for i = 1:4
    X = X_data{i}; Y = Y_data{i};
    
    % 使用 partialcorr 控制性别
    [r_part_sex, p_part_sex] = partialcorr(X, Y, FA_sex);
    delta_r = r_part_sex - base_r(i);
    
    fprintf('%-10s: partial r = %+.3f, p = %.3f (Δr = %+.3f)\n', vars{i}, r_part_sex, p_part_sex, delta_r);
end
fprintf('========================================================\n');