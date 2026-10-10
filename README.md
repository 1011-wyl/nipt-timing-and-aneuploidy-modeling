# NIPT Timing and Aneuploidy Modeling

基于 MATLAB 的 NIPT（无创产前检测）数学建模项目，围绕 Y 染色体浓度、孕妇 BMI、最佳检测时点和胎儿染色体异常判别展开分析。

## 项目内容

| 模块 | 主要任务 | 方法 |
| --- | --- | --- |
| 第一问 | 分析 Y 染色体浓度与孕周、BMI 的关系 | 相关分析、线性混合效应模型、残差分析 |
| 第二问 | 根据 BMI 分组并确定建议 NIPT 时点 | K-means、Kaplan–Meier 生存分析、蒙特卡洛误差分析 |
| 第三问 | 综合年龄、身高、体重等因素优化检测时点 | Davies–Bouldin 指数、多元回归、Bootstrap |
| 第四问 | 判别胎儿染色体是否异常 | 逻辑回归、混淆矩阵、ROC/AUC、特征重要性分析 |

## 目录结构

```text
├─ 第一问/
│  ├─ code1.m
│  └─ d1.xlsx
├─ 第二问/
│  ├─ code2_2.m
│  └─ d2.xlsx
├─ 第三问/
│  ├─ code3_1.m
│  └─ d3.xlsx
└─ 第四问/
   ├─ code4_2.m
   └─ d4.xlsx
```

## 运行环境

- MATLAB
- Statistics and Machine Learning Toolbox

运行某一问时，请先将 MATLAB 当前文件夹切换到对应子目录，再执行其中的 `.m` 文件，以便程序正确读取同目录下的 Excel 数据。

## 说明

本仓库用于数学建模学习与方法复现，不构成任何医学诊断或临床建议。
