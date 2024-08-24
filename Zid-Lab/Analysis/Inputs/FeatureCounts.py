##input: argv[1] is output file argv[2] is input file argv[3] is chrom file
##output: 8 Col file, should be the input for polairty.py

import sys
import os
import pandas as pd
import numpy as np

print(sys.argv[3])
ext = 'txt'
chromfiles = []
for root, dirs, files in os.walk(sys.argv[3]):
    for file in files:
        if file.split('.')[-1] == ext:
            chromfiles.append(sys.argv[3] + file)
chromfiles_s = sorted(chromfiles, key=lambda x: int(x.split('D')[1].split('.')[0]))

inputf = sys.argv[2]
outputf = str(sys.argv[1]) + str(inputf.split('/')[-1])
df = pd.read_csv(inputf, sep='\t', header=None, names=[str(x) for x in range(14)])
df = df.iloc[:, 0:4]
df = df.dropna() #removing rows containing null values.
roman_dict = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16, 'MT': 17}

chroms = []
posits = []
counts = []
for i in range(len(df)):
    fragL = df.iloc[i, 0][-2:] #last two char of the name of the reads
    orien = df.iloc[i, 1] #the orientation
    posit = df.iloc[i, 3] #for 5' read, position.
    chromr = df.iloc[i, 2]
    chromn = roman_dict.get(chromr)
    if orien == '-':
            posit_new = int(posit) + int(fragL)
    else:
            posit_new = int(posit) + 1
    counts.append(1)
    chroms.append(chromn)
    posits.append(posit_new)

df['2'] = chroms
df['3'] = posits
df['0'] = counts
df = df.groupby(['1','2','3']).agg({'0':'sum'})
df = df.reset_index()

names = []
chroms = []
posits = []
rposits = []
counts = []
starts = []
ends = []
lengths = []
oriens = []
prev_l = 0
print("before for loop")
for chromfile in chromfiles_s:
    chrom_num = chromfile.split('D')[1].split('.')[0]
    
    df_genes = pd.read_csv(chromfile, sep='\t', header=None, names=[str(x) for x in range(5)])
    df_genes = df_genes.set_index('0')
    #print(df_genes)
    #print(df_genes.columns)
    #print("in for loop")
    #df_genes are the chrom files, df is the input file.
    df_genes['3'] = pd.to_numeric(df_genes['3']) #df_genes['3'] is the end position of a gene
    df_genes['2'] = pd.to_numeric(df_genes['2']) #df_genes['3'] is the start position of a gene
    df_chrom = df[df['2'] == int(chrom_num)] #filter df to include only the rows that df['2'] == int(chrom_num) is true.
    print('Chromosome ',chrom_num)
    for i in range(len(df_chrom)):
        orien = df_chrom.iloc[i, -4]
        posit_new = df_chrom.iloc[i, -2]
        count = df_chrom.iloc[i, -1]
        #print(orien,posit_new,count)
        if orien == '-':                              # for - strand, the end posit should be greater than start posit -16 && the start posit should be less than 
            name_c = df_genes[(df_genes['4'] == orien) & (df_genes['3'] > (posit_new - 16)) & (df_genes['2'] < (posit_new - 13))].index.tolist()
            if len(name_c) > 0:
                name = name_c[0] 
                start = df_genes.loc[name, '3'] + 16
                end = df_genes.loc[name, '2'] + 13
                rposit = start - posit_new
        else:
            name_c = df_genes[(df_genes['4'] == orien) & (df_genes['3'] > (posit_new + 13)) & (df_genes['2'] < (posit_new + 16))].index.tolist()
            if len(name_c) > 0:
                name = name_c[0]
                start = df_genes.loc[name, '2'] - 16
                end = df_genes.loc[name, '3'] - 13
                rposit = posit_new - start #relative position is the relative position on the chromosome? 
        if len(name_c) > 0:
            names.append(name)
            starts.append(start)
            ends.append(end)
            lengths.append(abs(start - end)) 
            chroms.append(chrom_num)
            posits.append(posit_new)
            rposits.append(rposit)
            counts.append(count)
            oriens.append(orien)
    print('Total reads: ', str(len(df_chrom)))
    print('Aligned to genes: ', str(len(names)-prev_l))
    prev_l = len(names)


d = {'Name': names, 'Orient': oriens, 'Start': starts, 'End': ends, 'Position': posits, 'Length': lengths,'RelativePosition': rposits, 'Counts': counts}
dfcsv = pd.DataFrame(data=d)
dfcsv = dfcsv.sort_values(by=['Name'])
dfcsv = dfcsv.reset_index(drop=True)
dfcsv.to_csv(outputf)
