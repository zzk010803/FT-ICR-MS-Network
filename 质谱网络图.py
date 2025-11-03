#质谱网络图
#节点为分子式，边为可能发生的反应
# 首先导入我们需要的库
import pandas as pd
import numpy as np
import numpy.linalg as la
import pykrev as pk
import random
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# 读取示例数据集
dataset = pd.read_csv('D:\data\PBAT.csv')
msTuple = pk.read_corems(dataset)
intensityArray = msTuple.intensity
intensityNorm = pk.normalise_intensity(intensityArray)
formulaList = msTuple.formula
formulaMass = pk.calculate_mass(formulaList) # 计算数据集中每个公式的精确单同位素质量
#定义可能发生的反应
N = len(formulaList)
"""The reactions that can occur are defined below """
reactionDict = {
## 脱羧
'decarboxylation': -pk.calculate_mass(['CO2']),
## 甲基化
'methylation': pk.calculate_mass(['CH2']),
## 去甲基化
'demethylation': -pk.calculate_mass(['CH2']),
## 氢化作用
'hydrogenation': pk.calculate_mass(['H2']),
## 脱氢
'dehydrogenation': -pk.calculate_mass(['H2']),
## 水合 
'hydration': pk.calculate_mass(['H2O']),
## 脱水
'dehydration': -pk.calculate_mass(['H2O']),
# 氧化 
'oxidation': pk.calculate_mass(['O']),
# 还原
'reduction': -pk.calculate_mass(['O'])
}
r = pk.page_rank(msTuple, reactionDict = reactionDict, d = 0.9, tol = 0.01)
## 创建矩阵 L
L = np.zeros([N,N])
for j in range(N):
    for reactionType in reactionDict.keys():
        # 查找可能的匹配项，使用 np.round 来解释舍入误差
        matchIndex = np.where(np.round(formulaMass,8) == np.round((formulaMass[j] + reactionDict[reactionType]),8))
        # 如果找到匹配项
        if len(matchIndex[0]) > 0:
            L[matchIndex,j] = 1
    ## 对概率进行标准化，使它们之和为 1
    if sum(L[:,j]) > 0:
        L[:,j] = L[:,j]/sum(L[:,j])
    ## 如果列总和为零，则设置为 1/N
    else: 
        L[:,j] = 1/N
#运行 pageRank 算法以收敛
r = 100 * np.ones(N) / N # 设置概率向量
d = .9 # 阻尼系数 - 运行一次代码后，可以随意使用此参数。
M = d * L + (1-d)/N * np.ones([N, N]) # np.ones（） 是 J 矩阵，每个条目都有 1。
lastR = r
r = M @ r
i = 0
while la.norm(lastR - r) > 0.01 :
    lastR = r
    r = M @ r
    i += 1
print(str(i) + " iterations to convergence.")
argmax = np.argmax(r)
print('The highest page rank formula is')
print(formulaList[argmax])
print('With a probability at convergence of')
print(max(r))
print('The top 10 page rank formula are: ')
top_10_idx = np.argsort(r)[-10:]
for i in top_10_idx:
    print(formulaList[i])
print(f'The correlation between pagerank and peak intensity is')
print(stats.spearmanr(r,intensityNorm)[0])
print(f'The correlation between pagerank and formula mass is')
print(stats.spearmanr(r,formulaMass)[0])
print(f'A histogram of probabilities at d == {d}')
hist = plt.hist(r, bins = 100)
G, reactionCounts = pk.reaction_network(msTuple, 
                                        reactionDict = reactionDict, 
                                        filePath = 'PBAT1.gexf', 
                                        fileFormat = 'gexf', 
                                        nodeAnnotations = {"pageRank" : r, 
                                                           "peakIntensity" : intensityArray},
                                        roundVal = 8)
