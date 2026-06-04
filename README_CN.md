# FDRS：数字随机性筛查框架

FDRS（Fabrication-risk Digit Randomness Screening）是一个用于筛查原始数值型科研数据中非随机小数位数字模式的统计与机器学习辅助框架。

**注意：FDRS 只能作为辅助筛查和进一步核查优先级排序工具，不能作为判断数据伪造、篡改或科研不端的独立证据。**

## 主要功能

- 单个小数位数字分布检验；
- 联合两位小数位组合分布检验；
- 卡方统计量、P 值、Cramér’s V、Shannon entropy、标准化熵、KL divergence、标准化残差；
- 数字偏好指数；
- 渐进式抽样稳定性分析；
- 半监督机器学习风险评分；
- Random Forest、Elastic-net Logistic Regression、SVM radial、Isolation Forest 和 Ensemble 模型；
- 红蓝主色调可视化输出。

## 输入格式

输入文件为单列 txt 数值文件，每行一个独立数值：

```text
0.3167
0.4281
0.5029
0.3915
```

## 使用方式

先安装依赖：

```r
source("install_dependencies.R")
```

然后修改脚本开头的工作目录和输入文件名：

```r
setwd("YOUR/LOCAL/PROJECT/PATH")
INPUT_FILES <- c("RawData.txt", "ErrData.txt")
```

或：

```r
TARGET_FILES <- c("RealRawData2.txt")
```

再运行对应 R 脚本即可。

## 解释边界

高风险评分仅提示该数据集存在数字结构异常，需要进一步核查原始仪器导出文件、实验记录、数据处理流程、舍入规则、重复实验和专家审查。不能将 FDRS 输出直接等同于“造假”或“科研不端”。
