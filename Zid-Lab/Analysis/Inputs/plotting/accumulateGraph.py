# argv[1] input, WT afterFeatureCount file
# argv[2] input, MUTANT afterFeatureCount file
# argv[3] name of the gene
# argv[4] afterpolarity input file

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

name = sys.argv[3] #name of the gene
countf = sys.argv[4]
countdf = pd.read_csv(countf, sep=',', header=0) 
totalCount = 0
for row in range(len(countdf)):
    if countdf.loc[row, 'Name'] == name:
        totalCount = countdf.loc[row, 'Counts']
        break


inputf1 = sys.argv[1]
df1 = pd.read_csv(inputf1, sep=',', header=0)
for row in range(len(df1)):
    if df1.loc[row, "Name"] == name:
        length1 = df1.loc[row, "Length"]
        break
print(length1)
rposit = [0] * (length1 + 1)
for element in range(len(rposit)):
    rposit[element] = element
count1 = [0] * (length1 + 1)

for row in range(len(df1)):
    if df1.loc[row, 'Name'] == name:
        count1[df1.loc[row, 'RelativePosition']] = df1.loc[row, 'Counts'] / totalCount
for row in range(1, len(count1)):
    count1[row] = count1[row] + count1 [row - 1]


inputf2 = sys.argv[2]
df2 = pd.read_csv(inputf2, sep=',', header=0)
count2 = [0] * (length1 + 1)
for row in range(len(df2)):
    if df2.loc[row, 'Name'] == name:
        count2[df2.loc[row, 'RelativePosition']] = df2.loc[row, 'Counts'] / totalCount
        
        
for row in range(1, len(count1)):
    count2[row] = count2[row] + count2[row - 1]

WTname = os.path.splitext(os.path.basename(sys.argv[1]))[0]
mutantname = os.path.splitext(os.path.basename(sys.argv[2]))[0]
plt.ylabel("Counts", fontsize = 11)
plt.xlabel("Position", fontsize = 11)
plt.title(WTname + mutantname + name, fontsize = 13)
plt.plot(rposit, count1, label="WT")
plt.plot(rposit, count2, label="eIF5A Dep")
plt.legend(loc="lower right")
plt.show()