# input is argv[1] (the output file of featureCounts) 
# output is argv[2] (a 2 column file contains relativeposition and counts).

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

inputf = sys.argv[1]
outputf = sys.argv[2]
df = pd.read_csv(inputf, sep=',', header=0)
rposit = []
count = []
for row in range(len(df)):
    if df.loc[row, "Name"] == 'YPL063W':
        rposit.append(df.loc[row, "RelativePosition"])
        count.append(df.loc[row, "Counts"])
    

d = {"RelativePosition":rposit, "Counts":count}
output = pd.DataFrame(data=d)
output.to_csv(outputf)

plt.ylabel("Counts", fontsize = 11)
plt.xlabel("Position", fontsize = 11)
plt.title("YPL063W (137, dep)", fontsize = 13)

plt.bar(rposit, count, color="b", width=2.5)
plt.show()



