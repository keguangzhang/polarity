import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input1 = sys.argv[1] #T30
input2 = sys.argv[2] #Degron T0
input3 = sys.argv[3] #D208A
input4 = sys.argv[4] #WT
outputf = sys.argv[5] #output aggregated file

print('T30')
df1 = pd.read_csv(input1, sep=',', header=0)
T30Name = []
T30Length = []
T30Counts = []
T30Polarity = []
T30Space = []
for row in range(len(df1)):
    T30Name.append(df1.loc[row, 'Name'])
    T30Length.append(df1.loc[row, 'Length'])
    T30Counts.append(df1.loc[row, 'Counts'])
    T30Polarity.append(df1.loc[row, 'PolarityPerGene'])
    T30Space.append("\t")

df2 = pd.read_csv(input2, sep=',', header=0)
T0Name = ["None"] * len(df1)
T0Length = ["None"] * len(df1)
T0Counts = ["None"] * len(df1)
T0Polarity = ["None"] * len(df1)
T0Space = ["\t"] * len(df1)
notMatched1 = [] #records indexes of gene not matched to T30.
print('T0')
for row in range(len(df2)):
    if df2.loc[row, 'Name'] in df1.values:
        ind = df1.index[df1['Name']==df2.loc[row, 'Name']].tolist()
        T0Name[ind[0]] = df2.loc[row, 'Name']
        T0Length[ind[0]] = df2.loc[row, 'Length']
        T0Counts[ind[0]] = df2.loc[row, 'Counts']
        T0Polarity[ind[0]] = df2.loc[row, 'PolarityPerGene']
        
    else:
        notMatched1.append(row)
print("T0 not matched: " + str(len(notMatched1)) + '\n')

print('D208')
df3 = pd.read_csv(input3, sep=',', header=0)
D208Name = ["None"] * len(df1)
D208Length = ["None"] * len(df1)
D208Counts = ["None"] * len(df1)
D208Polarity = ["None"] * len(df1)
D208Space = ["\t"] * len(df1)
notMatched2 = []
for row in range(len(df3)):
    if df3.loc[row, 'Name'] in df1.values:
        ind = df1.index[df1['Name']==df3.loc[row, 'Name']].tolist()
        D208Name[ind[0]] = df3.loc[row, 'Name']
        D208Length[ind[0]] = df3.loc[row, 'Length']
        D208Counts[ind[0]] = df3.loc[row, 'Counts']
        D208Polarity[ind[0]] = df3.loc[row, 'PolarityPerGene']
        
    else:
        notMatched2.append(row)
print("D208 not matched: " + str(len(notMatched2)) + '\n')

print('WT')
df4 = pd.read_csv(input4, sep=',', header=0)
WTName = ["None"] * len(df1)
WTLength = ["None"] * len(df1)
WTCounts = ["None"] * len(df1)
WTPolarity = ["None"] * len(df1)
WTSpace = ["\t"] * len(df1)
notMatched3 = []
for row in range(len(df4)):
    if df4.loc[row, 'Name'] in df1.values:
        ind = df1.index[df1['Name']==df4.loc[row, 'Name']].tolist()
        WTName[ind[0]] = df4.loc[row, 'Name']
        WTLength[ind[0]] = df4.loc[row, 'Length']
        WTCounts[ind[0]] = df4.loc[row, 'Counts']
        WTPolarity[ind[0]] = df4.loc[row, 'PolarityPerGene']
        
    else:
        notMatched3.append(row)
print("WT not matched: " + str(len(notMatched3)) + '\n')

unmatchedName = []

uT30Length = ['None'] * 88
uT30Count = ['None'] * 88
uT30Polarity = ['None'] * 88
uT30Space = ['\t'] * 88
T30Length = T30Length + uT30Length
T30Counts = T30Counts + uT30Count
T30Polarity = T30Polarity + uT30Polarity 
T30Space = T30Space + uT30Space

uT0Length = ['None'] * 88
uT0Count = ['None'] * 88
uT0Polarity = ['None'] * 88
uT0Space = ['\t'] * 88
count = 0
for num in notMatched1:
    unmatchedName.append(df2.loc[num, 'Name'])
    uT0Length[count] = df2.loc[num, 'Length']
    uT0Count[count] = df2.loc[num, 'Counts']
    uT0Polarity[count] = df2.loc[num, 'PolarityPerGene']
    count = count + 1
T0Length = T0Length + uT0Length
T0Counts = T0Counts + uT0Count
T0Polarity = T0Polarity + uT0Polarity
T0Space = T0Space + uT0Space

uD208Length = ['None'] * 88
uD208Count = ['None'] * 88
uD208Polarity = ['None'] * 88  
uD208Space = ['\t'] * 88
for num in notMatched2:
    if df3.loc[num, 'Name'] not in unmatchedName:
        unmatchedName.append(df3.loc[num, 'Name'])
        uD208Length[count] = df3.loc[num, 'Length']
        uD208Count[count] = df3.loc[num, 'Counts']
        uD208Polarity[count] = df3.loc[num, 'PolarityPerGene']
        count = count + 1
    else:
        ind = unmatchedName.index(df3.loc[num, 'Name'])
        uD208Length[ind] = df3.loc[num, 'Length']
        uD208Count[ind] = df3.loc[num, 'Counts']
        uD208Polarity[ind] = df3.loc[num, 'PolarityPerGene']
D208Length = D208Length + uD208Length
D208Counts = D208Counts + uD208Count
D208Polarity = D208Polarity + uD208Polarity
D208Space = D208Space + uD208Space

uWTLength = ['None'] * 88
uWTCount = ['None'] * 88
uWTPolarity = ['None'] * 88
uWTSpace = ['\t'] * 88
for num in notMatched3:
    if df4.loc[num, 'Name'] not in unmatchedName:
        unmatchedName.append(df4.loc[num, 'Name'])
        uWTLength[count] = df4.loc[num, 'Length']
        uWTCount[count] = df4.loc[num, 'Counts']
        uWTPolarity[count] = df4.loc[num, 'PolarityPerGene']
        count = count + 1
    else:
        ind = unmatchedName.index(df4.loc[num, 'Name'])
        uWTLength[ind] = df4.loc[num, 'Length']
        uWTCount[ind] = df4.loc[num, 'Counts']
        uWTPolarity[ind] = df4.loc[num, 'PolarityPerGene']
WTLength = WTLength + uWTLength
WTCounts = WTCounts + uWTCount
WTPolarity = WTPolarity + uWTPolarity
WTSpace = WTSpace + uWTSpace

T30Name = T30Name + unmatchedName
print('T30: '+str(len(T30Name))+'\t'+str(len(T30Length))+ '\t'+str(len(T30Counts))+'\t'+str(len(T30Polarity))+'\t'+str(len(T30Space))+'\n')
print('T0: '+str(len(T0Length))+ '\t'+str(len(T0Counts))+'\t'+str(len(T0Polarity))+'\t'+str(len(T0Space))+'\n')
print('D208: '+'\t'+str(len(D208Length))+ '\t'+str(len(D208Counts))+'\t'+str(len(D208Polarity))+'\t'+str(len(D208Space))+'\n')
print('WT: '+'\t'+str(len(WTLength))+ '\t'+str(len(WTCounts))+'\t'+str(len(WTPolarity))+'\n')
print(len(unmatchedName))
print(count)


d = {"Name":T30Name, "D208ALength":D208Length, "D208ACounts":D208Counts, "D208APolarity":D208Polarity, "\t":D208Space,
     "T0Length":T0Length, "T0Counts":T0Counts, "T0Polarity":T0Polarity, "\t\t":T0Space,
     "T30Length":T30Length, "T30Counts":T30Counts, "PT30olarity":T30Polarity, "\t\t\t":T30Space,
     "WTLength":WTLength, "WTCounts":WTCounts, "WTPolarity":WTPolarity}
output = pd.DataFrame(data=d)
output.sort_values('Name')
output.to_csv(outputf)








