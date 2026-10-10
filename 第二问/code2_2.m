data = readtable('d2.xlsx'); 

% 简单数据处理
uniqueIDs = unique(data.MotherID);
nMothers = length(uniqueIDs);
eventTime = zeros(nMothers, 1);
eventStatus = zeros(nMothers, 1); % 1: 达标, 0: 未达标
bmiValues = zeros(nMothers, 1);

% 计算每位孕妇的首次达标时间
for i = 1:nMothers
    id = uniqueIDs{i};
    subData = data(strcmp(data.MotherID, id), :);
    [sortedGA, idx] = sort(subData.GestationalAge);
    sortedYconc = subData.Yconc(idx);
    
    % 找到首次达标时间
    idx_event = find(sortedYconc >= 0.04, 1);
    if ~isempty(idx_event)
        eventTime(i) = sortedGA(idx_event);
        eventStatus(i) = 1;
    else
        eventTime(i) = sortedGA(end);
        eventStatus(i) = 0;
    end
    bmiValues(i) = subData.BMI(1);
end

% 由拐点可得3为最佳分组数
optimal_k = 3;

% 使用最佳K值进行K-means聚类
[idx_bmi, centroids] = kmeans(bmiValues, optimal_k);

% 确定每组的BMI范围
groups = unique(idx_bmi);
nGroups = length(groups);
group_bmi_range = cell(nGroups, 1);

for i = 1:nGroups
    bmi_i = bmiValues(idx_bmi == groups(i));
    group_bmi_range{i} = [min(bmi_i), max(bmi_i)];
    fprintf('组%d: BMI范围 [%.2f, %.2f], 数量%d\n', i, group_bmi_range{i}(1), group_bmi_range{i}(2), length(bmi_i));
end

% 对每个BMI分组进行生存分析，计算最佳NIPT时点
group_t95 = zeros(nGroups, 1);

for i = 1:nGroups
    group_i = (idx_bmi == groups(i));
    time_i = eventTime(group_i);
    status_i = eventStatus(group_i);
    
    % 使用Kaplan-Meier方法估计生存函数
    [kmf, x] = ecdf(time_i, 'Censoring', 1-status_i, 'Function', 'survivor');
    
    % 找到生存函数值首次下降到0.05以下的时间点
    idx_t95 = find(kmf <= 0.05, 1);
    if isempty(idx_t95)
        group_t95(i) = x(end);
    else
        group_t95(i) = x(idx_t95);
    end
    
    fprintf('组%d: 最佳NIPT时点: %.2f周\n', i, group_t95(i));
end

% 误差分析
nSimulations = 100; % 模拟次数
sigma_error = 0.005; % 误差标准差
group_t95_error = zeros(nGroups, nSimulations);

for sim = 1:nSimulations
    eventTime_error = zeros(nMothers, 1);
    eventStatus_error = zeros(nMothers, 1);
    
    for i = 1:nMothers
        id = uniqueIDs{i};
        subData = data(strcmp(data.MotherID, id), :);
        [sortedGA, idx] = sort(subData.GestationalAge);
        
        % 添加随机误差到Y浓度
        sortedYconc = subData.Yconc(idx) + normrnd(0, sigma_error, size(subData.Yconc(idx)));
        
        % 找到首次达标时间
        idx_event = find(sortedYconc >= 0.04, 1);
        if ~isempty(idx_event)
            eventTime_error(i) = sortedGA(idx_event);
            eventStatus_error(i) = 1;
        else
            eventTime_error(i) = sortedGA(end);
            eventStatus_error(i) = 0;
        end
    end
    
    % 计算每个分组的最佳NIPT时点
    for i = 1:nGroups
        group_i = (idx_bmi == groups(i));
        time_i = eventTime_error(group_i);
        status_i = eventStatus_error(group_i);
        
        [kmf, x] = ecdf(time_i, 'Censoring', 1-status_i, 'Function', 'survivor');
        
        idx_t95 = find(kmf <= 0.05, 1);
        if isempty(idx_t95)
            group_t95_error(i, sim) = x(end);
        else
            group_t95_error(i, sim) = x(idx_t95);
        end
    end
end

% 计算误差情况下最佳NIPT时点的均值和95%置信区间
group_t95_mean = mean(group_t95_error, 2);
group_t95_ci = zeros(nGroups, 2);
for i = 1:nGroups
    sorted_t95 = sort(group_t95_error(i, :));
    group_t95_ci(i, 1) = sorted_t95(round(0.025 * nSimulations));
    group_t95_ci(i, 2) = sorted_t95(round(0.975 * nSimulations));
end

disp('误差分析结果:');
for i = 1:nGroups
    fprintf('组%d: 无误差时点: %.2f周, 有误差时点均值: %.2f周, 95%%置信区间: [%.2f, %.2f]周\n', ...
        i, group_t95(i), group_t95_mean(i), group_t95_ci(i, 1), group_t95_ci(i, 2));
end

figure;
hold on;
for i = 1:nGroups
    errorbar(i, group_t95_mean(i), group_t95_mean(i) - group_t95_ci(i, 1), ...
        group_t95_ci(i, 2) - group_t95_mean(i), 'o', 'LineWidth', 2);
    plot(i, group_t95(i), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end
xlabel('BMI分组');
ylabel('最佳NIPT时点(周)');
title('误差对最佳NIPT时点的影响');
legend('有误差(均值±95% CI)', '无误差');
set(gca, 'XTick', 1:nGroups);
grid on;
hold off;