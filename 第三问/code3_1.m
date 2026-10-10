% 首先对BMI进行聚类，然后考虑其他因素
clear; clc; close all;

data = readtable('d3.xlsx');
rng(99);
 
% 数据预处理
motherIDs = unique(data.MotherID);
numMothers = length(motherIDs);

% 初始化变量
bmi = zeros(numMothers, 1);
ages = zeros(numMothers, 1);
heights = zeros(numMothers, 1);
weights = zeros(numMothers, 1);
eventTimes = zeros(numMothers, 1);

% 计算检测误差
errorVariances = [];
for i = 1:numMothers
    motherData = data(strcmp(data.MotherID, motherIDs{i}), :);
    
    % 提取每位母亲的特征
    bmi(i) = motherData.BMI(1);
    ages(i) = motherData.Age(1);
    heights(i) = motherData.Height(1);
    weights(i) = motherData.Weight(1);
    
    % 检查是否有重复测量
    [uniqueGA, ~, ic] = unique(motherData.GestationalAge);
    for j = 1:length(uniqueGA)
        sameGA = motherData(ic == j, :);
        if height(sameGA) > 1
            % 计算同一孕周测量的方差
            errorVariances = [errorVariances; var(sameGA.Yconc)];
        end
    end
    
    % 估计达标时间
    if any(motherData.Yconc >= 0.04)
        % 如果观测到达标，记录首次达标的时间
        firstEventIndex = find(motherData.Yconc >= 0.04, 1);
        eventTimes(i) = motherData.GestationalAge(firstEventIndex);
    else
        % 如果未观测到达标，使用线性拟合估计达标时间
        if length(motherData.GestationalAge) > 1
            % 确保有足够的数据点进行拟合
            p = polyfit(motherData.GestationalAge, motherData.Yconc, 1);
            if p(1) > 0 % 只有正斜率才有意义
                eventTimes(i) = (0.04 - p(2)) / p(1);
            end
        end
    end
end

% 计算平均检测误差方差
avgErrorVariance = mean(errorVariances);
sigma = sqrt(avgErrorVariance);
fprintf('估计的检测误差标准差: %.6f\n', sigma);

% 移除无法估计达标时间的样本
validIndices = ~isnan(eventTimes);
bmi = bmi(validIndices);
ages = ages(validIndices);
heights = heights(validIndices);
weights = weights(validIndices);
eventTimes = eventTimes(validIndices);
motherIDs = motherIDs(validIndices);
numMothers = sum(validIndices);

fprintf('有效样本数量: %d\n', numMothers);

% 计算Davies-Bouldin指数
kList = 2:4;
dbIndex = zeros(length(kList), 1);

figure;
for i = 1:length(kList)
    k = kList(i);
    [idx, centers] = kmeans(bmi, k);
    dbIndex(i) = db_index(bmi, idx, centers);
    % 绘制聚类结果
    subplot(1, length(kList), i);
    gscatter(1:numMothers, bmi, idx);
    title(sprintf('K = %d, DB指数 = %.4f', k, dbIndex(i)));
    xlabel('样本索引');
    ylabel('BMI');
end

% 选择最佳聚类数量
[~, bestKIdx] = min(dbIndex);
bestK = kList(bestKIdx);
fprintf('最佳聚类数量: %d\n', bestK);

% 使用最佳K值进行最终聚类
[idx, centers] = kmeans(bmi, bestK);

% 显示聚类结果
fprintf('\nBMI聚类结果:\n');
for i = 1:bestK
    groupBmi = bmi(idx == i);
    groupMothers = motherIDs(idx == i);
    
    fprintf('组%d (平均BMI=%.2f, 范围=[%.1f, %.1f],数量%d): ', i, ...
        mean(groupBmi), min(groupBmi), max(groupBmi), length(groupBmi));
    for j = 1:length(groupMothers)
        fprintf('%s ', groupMothers{j});
    end
    fprintf('\n');
end

% 对于每个BMI组，分析年龄、身高、体重的分布
fprintf('\n各组特征分析:\n');
for i = 1:bestK
    groupIndices = (idx == i);
    groupAges = ages(groupIndices);
    groupHeights = heights(groupIndices);
    groupWeights = weights(groupIndices);
    groupEventTimes = eventTimes(groupIndices);
    fprintf('组%d:\n', i);
    fprintf('  年龄: 平均值=%.1f, 范围=[%.1f, %.1f]\n', ...
            mean(groupAges), min(groupAges), max(groupAges));
    fprintf('  身高: 平均值=%.1f, 范围=[%.1f, %.1f]\n', ...
            mean(groupHeights), min(groupHeights), max(groupHeights));
    fprintf('  体重: 平均值=%.1f, 范围=[%.1f, %.1f]\n', ...
            mean(groupWeights), min(groupWeights), max(groupWeights));
    fprintf('  达标时间: 平均值=%.2f周, 范围=[%.2f, %.2f]周\n', ...
            mean(groupEventTimes), min(groupEventTimes), max(groupEventTimes));
end

% 对于每个BMI组，建立达标时间与年龄、身高、体重的回归模型
fprintf('\n多因素回归分析:\n');
for i = 1:bestK
    groupIndices = (idx == i);
    if sum(groupIndices) > 2 
        X = [ones(sum(groupIndices), 1), ages(groupIndices), heights(groupIndices), weights(groupIndices)];
        y = eventTimes(groupIndices);
        [b, bint, r, rint, stats] = regress(y, X);
        fprintf('组%d回归模型:\n', i);
        fprintf('  截距: %.3f (95%% CI: [%.3f, %.3f])\n', b(1), bint(1,1), bint(1,2));
        fprintf('  年龄系数: %.3f (95%% CI: [%.3f, %.3f])\n', b(2), bint(2,1), bint(2,2));
        fprintf('  身高系数: %.3f (95%% CI: [%.3f, %.3f])\n', b(3), bint(3,1), bint(3,2));
        fprintf('  体重系数: %.3f (95%% CI: [%.3f, %.3f])\n', b(4), bint(4,1), bint(4,2));
        fprintf('  R² = %.3f, F统计量 = %.3f, p值 = %.3f\n', stats(1), stats(2), stats(3));
    end
end

% 计算每组的最佳NIPT时点
bestTimes = zeros(bestK, 1);
for i = 1:bestK
    groupTimes = eventTimes(idx == i);
    bestTimes(i) = prctile(groupTimes, 95);
    fprintf('组%d的最佳NIPT时点: %.2f周\n', i, bestTimes(i));
end

% 使用Bootstrap方法模拟检测误差对结果的影响
numBoot = 1000; % Bootstrap次数
bootBestTimes = zeros(numBoot, bestK); 

for b = 1:numBoot
    % 添加随机误差到Yconc
    dataBoot = data;
    dataBoot.Yconc = dataBoot.Yconc + randn(height(dataBoot), 1) * sigma;
    
    % 重新估计达标时间
    bootEventTimes = zeros(numMothers, 1);
    for i = 1:numMothers
        motherData = dataBoot(strcmp(dataBoot.MotherID, motherIDs{i}), :);
        
        if any(motherData.Yconc >= 0.04)
            firstEventIndex = find(motherData.Yconc >= 0.04, 1);
            bootEventTimes(i) = motherData.GestationalAge(firstEventIndex);
        else
            if length(motherData.GestationalAge) > 1
                p = polyfit(motherData.GestationalAge, motherData.Yconc, 1);
                if p(1) > 0
                    bootEventTimes(i) = (0.04 - p(2)) / p(1);
                else
                    bootEventTimes(i) = NaN;
                end
            else
                bootEventTimes(i) = NaN;
            end
        end
    end
    
    % 移除无效值
    validBoot = ~isnan(bootEventTimes);
    bootEventTimes = bootEventTimes(validBoot);
    bootBmi = bmi(validBoot);
    bootIdx = idx(validBoot);
    
    % 计算每组的最佳时点
    for i = 1:bestK
        groupTimes = bootEventTimes(bootIdx == i);
        if ~isempty(groupTimes)
            bootBestTimes(b, i) = prctile(groupTimes, 95);
        else
            bootBestTimes(b, i) = NaN;
        end
    end
end

% 计算置信区间
alpha = 0.05; % 95%置信区间
ci = zeros(bestK, 2);
for i = 1:bestK
    validBoot = ~isnan(bootBestTimes(:, i));
    if sum(validBoot) > 0
        ci(i, :) = prctile(bootBestTimes(validBoot, i), [100*alpha/2, 100*(1-alpha/2)]);
        fprintf('组%d最佳时点的95%%置信区间: [%.2f, %.2f]周\n', i, ci(i, 1), ci(i, 2));
    end
end

% 绘制Bootstrap结果
figure;
for i = 1:bestK
    subplot(1, bestK, i);
    histogram(bootBestTimes(:, i), 20);
    title(sprintf('组%d最佳时点的Bootstrap分布', i));
    xlabel('最佳时点(周)');
    ylabel('频率');
    line([bestTimes(i), bestTimes(i)], ylim, 'Color', 'r', 'LineWidth', 2);
    legend('Bootstrap分布', '原始估计');
end

% 绘制BMI与达标时间的关系，按组着色
figure;
scatter(bmi, eventTimes, 100, idx, 'filled');
xlabel('BMI');
ylabel('达标时间(周)');
title('BMI与达标时间的关系（按组着色）');
colorbar;

% 绘制其他因素与达标时间的关系
figure;
subplot(2, 2, 1);
scatter(ages, eventTimes, 100, idx, 'filled');
xlabel('年龄');
ylabel('达标时间(周)');
title('年龄与达标时间的关系');

subplot(2, 2, 2);
scatter(heights, eventTimes, 100, idx, 'filled');
xlabel('身高(cm)');
ylabel('达标时间(周)');
title('身高与达标时间的关系');

subplot(2, 2, 3);
scatter(weights, eventTimes, 100, idx, 'filled');
xlabel('体重(kg)');
ylabel('达标时间(周)');
title('体重与达标时间的关系');

subplot(2, 2, 4);
scatter(bmi, eventTimes, 100, idx, 'filled');
xlabel('BMI');
ylabel('达标时间(周)');
title('BMI与达标时间的关系');

% 计算Davies-Bouldin指数
function db = db_index(data, idx, centers)
    k = size(centers, 1);
    n = size(data, 1);
    
    % 计算每个簇的直径（平均距离）
    S = zeros(k, 1);
    for i = 1:k
        clusterData = data(idx == i);
        if length(clusterData) > 1
            S(i) = mean(sqrt((clusterData - centers(i)).^2));
        else
            S(i) = 0;
        end
    end
    
    % 计算簇间距离
    R = zeros(k, k);
    for i = 1:k
        for j = i+1:k
            Mij = sqrt((centers(i) - centers(j))^2);
            R(i, j) = (S(i) + S(j)) / Mij;
        end
    end
    
    % 计算每个簇的最大R值
    D = zeros(k, 1);
    for i = 1:k
        if k > 1
            D(i) = max(R(i, setdiff(1:k, i)));
        else
            D(i) = 0;
        end
    end
    
    db = mean(D);
end