data = readtable('d4.xlsx');

features = data(:, {'Age', 'BMI', 'TotalReads', 'GCContent', ...
                   'Zscore13', 'Zscore18', 'Zscore21', 'Xz', 'Xconc', ...
                   'Chr13GCContent', 'Chr18GCContent', 'Chr21GCContent', ...
                   'FilteredReadRate', 'GestationalAge'});

% 将染色体异常转换为二进制目标变量
target = zeros(height(data), 1);
for i = 1:height(data)
    if isempty(data.ChromosomalAneuploidy{i})
        target(i) = 0; % 正常
    else
        target(i) = 1; % 异常
    end
end

% 处理缺失值
feature_array = table2array(features);
nan_indices = any(isnan(feature_array), 2);
feature_array = feature_array(~nan_indices, :);
target = target(~nan_indices);

% 标准化特征
feature_means = mean(feature_array);
feature_stds = std(feature_array);
normalized_features = (feature_array - feature_means) ./ feature_stds;

% 划分训练集和测试集
rng(1);
cv = cvpartition(size(normalized_features, 1), 'HoldOut', 0.3);
idx_train = training(cv);
idx_test = test(cv);

X_train = normalized_features(idx_train, :);
y_train = target(idx_train);
X_test = normalized_features(idx_test, :);
y_test = target(idx_test);

% 使用逻辑回归构建模型
mdl = fitglm(X_train, y_train, 'Distribution', 'binomial', 'Link', 'logit');
fprintf('逻辑回归模型详细信息\n');

fprintf('模型系数:\n');
feature_names = features.Properties.VariableNames;
fprintf('截距 (β₀): %.6f\n', mdl.Coefficients.Estimate(1));
for i = 1:length(feature_names)
    fprintf('%s (β%d): %.6f\n', feature_names{i}, i, mdl.Coefficients.Estimate(i+1));
end
fprintf('\n');

fprintf('异常概率计算公式:\n');
fprintf('P(异常) = 1 / (1 + exp(-(%.6f', mdl.Coefficients.Estimate(1));
for i = 1:length(feature_names)
    if mdl.Coefficients.Estimate(i+1) >= 0
        fprintf(' + %.6f*%s', mdl.Coefficients.Estimate(i+1), feature_names{i});
    else
        fprintf(' - %.6f*%s', abs(mdl.Coefficients.Estimate(i+1)), feature_names{i});
    end
end
fprintf(')))\n\n');

% 系数的统计显著性
fprintf('系数统计显著性 (p值):\n');
fprintf('截距: p = %.6f\n', mdl.Coefficients.pValue(1));
for i = 1:length(feature_names)
    fprintf('%s: p = %.6f', feature_names{i}, mdl.Coefficients.pValue(i+1));
    if mdl.Coefficients.pValue(i+1) < 0.05
        fprintf(' *'); 
    end
    fprintf('\n');
end
fprintf('* 表示p < 0.05，统计显著\n\n');

% 模型拟合优度
fprintf('模型拟合优度:\n');
fprintf('对数似然值: %.4f\n', mdl.LogLikelihood);
fprintf('AIC: %.4f\n', mdl.ModelCriterion.AIC);
fprintf('BIC: %.4f\n', mdl.ModelCriterion.BIC);
fprintf('偏差: %.4f\n', mdl.Deviance);

y_pred_prob = predict(mdl, X_test);
y_pred = double(y_pred_prob > 0.5);

result_table = table();
result_table.Real_Category = y_test;
result_table.Prediction_probability  = y_pred_prob;
result_table.Forecast_categories = y_pred;
result_table.Correct = (y_test == y_pred);

% 计算并输出模型性能指标
accuracy = sum(y_pred == y_test) / length(y_test);
confusion_mat = confusionmat(y_test, y_pred);

if size(confusion_mat, 1) == 2
    TP = confusion_mat(2,2);
    TN = confusion_mat(1,1);
    FP = confusion_mat(1,2);
    FN = confusion_mat(2,1);
    
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    specificity = TN / (TN + FP);
    f1_score = 2 * (precision * recall) / (precision + recall);
    
    fprintf('\n模型性能指标:\n');
    fprintf('准确率: %.2f%%\n', accuracy*100);
    fprintf('精确率: %.2f%%\n', precision*100);
    fprintf('召回率: %.2f%%\n', recall*100);
    fprintf('特异性: %.2f%%\n', specificity*100);
    fprintf('F1分数: %.2f%%\n', f1_score*100);
end

% 特征重要性图表
coefficients = mdl.Coefficients.Estimate(2:end);
[~, idx] = sort(abs(coefficients), 'descend');

figure;
barh(coefficients(idx));
set(gca, 'YTickLabel', feature_names(idx));
title('特征重要性（逻辑回归系数）');
xlabel('系数绝对值');
grid on;

% 绘制ROC曲线
if length(unique(y_test)) > 1
    [X_roc, Y_roc, T, AUC] = perfcurve(y_test, y_pred_prob, 1);
    
    figure;
    plot(X_roc, Y_roc, 'LineWidth', 2);
    hold on;
    plot([0 1], [0 1], '--', 'Color', [0.5 0.5 0.5]);
    xlabel('假阳性率');
    ylabel('真阳性率');
    title(sprintf('ROC曲线 (AUC = %.3f)', AUC));
    legend(sprintf('逻辑回归 (AUC=%.3f)', AUC), '随机分类', 'Location', 'southeast');
    grid on;
    hold off;
end

% 创建一个可用于新样本预测的函数
fprintf('新样本预测函数\n');
fprintf('function is_abnormal = predict_abnormal(patient_data)\n');
fprintf('     标准化新数据\n');
fprintf('    normalized_data = (patient_data - [');
fprintf('%.6f ', feature_means);
fprintf(']) ./ [');
fprintf('%.6f ', feature_stds);
fprintf('];\n');
fprintf('    \n');
fprintf('     计算线性组合\n');
fprintf('    linear_combination = %.6f', mdl.Coefficients.Estimate(1));
for i = 1:length(feature_names)
    if mdl.Coefficients.Estimate(i+1) >= 0
        fprintf(' + %.6f * normalized_data(%d)', mdl.Coefficients.Estimate(i+1), i);
    else
        fprintf(' - %.6f * normalized_data(%d)', abs(mdl.Coefficients.Estimate(i+1)), i);
    end
end
fprintf(';\n');
fprintf('    \n');
fprintf('    计算异常概率\n');
fprintf('    prob = 1 / (1 + exp(-linear_combination));\n');
fprintf('    \n');
fprintf('    判定结果\n');
fprintf('    is_abnormal = prob > 0.5;\n');
fprintf('end\n');