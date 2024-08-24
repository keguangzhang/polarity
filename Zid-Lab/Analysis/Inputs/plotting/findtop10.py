#argv[1] input file gives the counts and relative positions on the gene.
#argv[2] input file gives the total number of ribosomes on the gene.
#(useless, just to see) argv[3] output file is the all the top ten gene counts and relative position.


import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input1 = sys.argv[1]
input2 = sys.argv[2]
outputf = sys.argv[3]

#find the 10 most abundunt gene. Put the gene counts in array top10num. Put the gene names in array top10name.
df = pd.read_csv(input2, sep=',', header=0)
# an array record the name of the top 10 genes
top10num = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top10name = [None, None, None, None, None, None, None, None, None, None]
for num in range(10): 
    for row in range(len(df)):
        if df.loc[row, "Counts"] > top10num[num] & df.loc[row, "Length"] > 1500:
            found = False
            for name in range(len(top10name)):
                if df.loc[row, "Name"] == top10name[name]:
                    found = True
            if found == False:
                top10num[num] = df.loc[row, "Counts"] 
                top10name[num] = df.loc[row, "Name"]


# read the output file of the output of the featurecount to find the count and relative positions of the top 10 genes.
dataframe = pd.read_csv(input1, sep=',', header=0)
names = []
rposit = []
count = []
for num in range(10):
    for row in range(len(dataframe)):
        if dataframe.loc[row, "Name"] == top10name[num]:
            names.append(top10name[num])
            rposit.append(dataframe.loc[row, "RelativePosition"])
            count.append(dataframe.loc[row, "Counts"])

d = {"Name": names, "RelativePosition":rposit, "Counts":count}
output = pd.DataFrame(data=d)
output.to_csv(outputf)


#Standardize counts: read each line in top10.csv. Calculate count/total count of the entire gene 
#and then change the count of each line into the calculated value.
df2 = pd.read_csv(sys.argv[3])
rcount = []
for row in range(len(df2)):
    for name in range(len(top10name)):
        if df2.loc[row, "Name"] == top10name[name]:
            rcount.append(float(df2.loc[row, "Counts"] / top10num[name]))

            
d2 = {"Name": names, "RelativePosition": rposit, "Count":count, "RelativeCount": rcount}
output = pd.DataFrame(data=d2)
output.to_csv(outputf) 

#Add counts of the same position together.
#Find the largest length

df3 = pd.read_csv(sys.argv[3])
col = "RelativePosition"
maxV = df3["RelativePosition"].max()

#create a new array with all 0 that has the length of the largest length
aggregated = []
for num2 in range(maxV+1):
    aggregated.append(0)


#add up all the counts with the same relative position
for row in range(len(df3)):
    aggregated[df3.loc[row, "RelativePosition"]] = aggregated[df3.loc[row, "RelativePosition"]] + df3.loc[row, "RelativeCount"]
    

rp = []
for num3 in range(maxV+1):
    rp.append(num3)

d3 = {"RelativePosition":rp, "TotalCount":aggregated}
output = pd.DataFrame(data=d3)
output.to_csv(outputf) 

plt.ylabel("Counts", fontsize = 11)
plt.xlabel("Position", fontsize = 11)
plt.title("top10(881)", fontsize = 13)

plt.bar(rp, aggregated, color="b", width=2.5)
plt.show()


