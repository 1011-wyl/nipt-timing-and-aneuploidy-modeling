clear; close all; clc;

data = readtable('d1.xlsx');

% 数据预处理 转换孕周格式
ga_numeric = zeros(height(data), 1);
for i = 1:height(data)
    ga_str = data.GestationalAge{i};
    % 使用正则表达式提取周数和天数
    tokens = regexp(ga_str, '(\d+)w\+(\d+)', 'tokens');
    if ~isempty(tokens)
        weeks = str2double(tokens{1}{1});
        days = str2double(tokens{1}{2});
        ga_numeric(i) = weeks + days/7;
    else
        tokens = regexp(ga_str, '(\d+)w', 'tokens');
        if ~isempty(tokens)
            weeks = str2double(tokens{1}{1});
            ga_numeric(i) = weeks;
        else
            ga_numeric(i) = NaN;
        end
    end
end

analysis_data = table(data.MotherID, ga_numeric, data.BMI, data.Yconc, ...
    'VariableNames', {'MotherID', 'GA', 'BMI', 'Yconc'});

analysis_data = rmmissing(analysis_data);

% 绘制散点图
figure;
subplot(1, 2, 1);
scatter(analysis_data.GA, analysis_data.Yconc, 'filled');
xlabel('孕周 (周)');
ylabel('Y染色体浓度');
title('Y染色体浓度与孕周的关系');
grid on;

subplot(1, 2, 2);
scatter(analysis_data.BMI, analysis_data.Yconc, 'filled');
xlabel('BMI');
ylabel('Y染色体浓度');
title('Y染色体浓度与BMI的关系');
grid on;

% 计算相关系数
corr_ga = corr(analysis_data.GA, analysis_data.Yconc);
corr_bmi = corr(analysis_data.BMI, analysis_data.Yconc);
fprintf('Y染色体浓度与孕周的相关系数: %.4f\n', corr_ga);
fprintf('Y染色体浓度与BMI的相关系数: %.4f\n', corr_bmi);

% 拟合线性混合效应模型
lme = fitlme(analysis_data, 'Yconc ~ GA + BMI + (1|MotherID)');

disp('线性混合效应模型结果:');
disp(lme);

% 残差分析
figure;
subplot(2, 2, 1);
plotResiduals(lme, 'fitted');
xlabel('拟合值');
ylabel('残差');
title('残差与拟合值的关系');

subplot(2, 2, 2);
plotResiduals(lme, 'probability');
title('残差正态概率图');

subplot(2, 2, 3);
plotResiduals(lme, 'lagged');
title('残差自相关图');

subplot(2, 2, 4);
histogram(residuals(lme));
xlabel('残差');
ylabel('频数');
title('残差分布');
y_pred = predict(lme);
figure;
plot(analysis_data.Yconc, y_pred, 'o');
hold on;
plot([min(analysis_data.Yconc), max(analysis_data.Yconc)], ...
     [min(analysis_data.Yconc), max(analysis_data.Yconc)], 'r--');
xlabel('实际Y染色体浓度');
ylabel('预测Y染色体浓度');
title('预测值与实际值的关系');
grid on;
legend('数据点', '理想拟合线', 'Location', 'best');

% 计算模型评估指标
mse = mean((analysis_data.Yconc - y_pred).^2);
rmse = sqrt(mse);
r2 = 1 - sum((analysis_data.Yconc - y_pred).^2) / sum((analysis_data.Yconc - mean(analysis_data.Yconc)).^2);
fprintf('模型评估指标:\n');
fprintf('均方误差 (MSE): %.6f\n', mse);
fprintf('均方根误差 (RMSE): %.6f\n', rmse);
fprintf('决定系数 (R²): %.4f\n', r2);

% 显著性检验
hypothesis_test = coefTest(lme);
fprintf('模型整体显著性检验 (p值): %.6f\n', hypothesis_test);

% 固定效应系数及其显著性
fixed_effects = lme.Coefficients;
disp('固定效应系数及其显著性:');
disp(fixed_effects);

% 随机效应方差
random_effects = lme.randomEffects;
fprintf('随机效应方差: %.6f\n', var(random_effects));