file1 = open('./alignedResults/finalResults137.fastq','r')
output = open('finalfinal', 'w')
AlignCountTotal = 0.0
for line in file1:
    temp = line.rstrip()
    temp1 = temp.split('\t') 
    Chr = int(temp1[1])
    Posit = int(temp1[3])
    Posit = Posit + 1
    Orien = temp1[0]
    count = temp1[2]
    AlignCount = 0
    
    #Seq = temp1[5]
    if Chr == 1 :
        file2 = open('./Chrom/ORF_ND1.txt','r')
    elif Chr == 2 :
        file2 = open('./Chrom/ORF_ND2.txt','r')
    elif Chr == 3 :
        file2 = open('./Chrom/ORF_ND3.txt','r')
    elif Chr == 4 :
        file2 = open('./Chrom/ORF_ND4.txt','r')
    elif Chr == 5 :
        file2 = open('./Chrom/ORF_ND5.txt','r')
    elif Chr == 6 :
        file2 = open('./Chrom/ORF_ND6.txt','r')
    elif Chr == 7 :
        file2 = open('./Chrom/ORF_ND7.txt','r')
    elif Chr == 8 :
        file2 = open('./Chrom/ORF_ND8.txt','r')
    elif Chr == 9 :
        file2 = open('./Chrom/ORF_ND9.txt','r')
    elif Chr == 10 :
        file2 = open('./Chrom/ORF_ND10.txt','r')
    elif Chr == 11 :
        file2 = open('./Chrom/ORF_ND11.txt','r')
    elif Chr == 12 :
        file2 = open('./Chrom/ORF_ND12.txt','r')
    elif Chr == 13 :
        file2 = open('./Chrom/ORF_ND13.txt','r')
    elif Chr == 14 :
        file2 = open('./Chrom/ORF_ND14.txt','r')
    elif Chr == 15 :
        file2 = open('./Chrom/ORF_ND15.txt','r')
    elif Chr == 16 :
        file2 = open('./Chrom/ORF_ND16.txt','r')
    elif Chr == 17 :
        file2 = open('./Chrom/ORF_ND17.txt','r')    
    
    for line2 in file2 :
        temp3 = line2.rstrip()
        temp2 = temp3.split('\t')
        Chr2 = int(temp2[1])
        if Chr == Chr2: #checking if the same chromosome
            if Orien == '-' :  # C orientation
                name = temp2[0]
                Orien2 = temp2[4]
                if Orien2 == '-' :  #also C orientation
                    tStrt2 = int(temp2[3])
                    tEnd2 = int(temp2[2])
                    Strt2 = tStrt2+16
                    End2 = tEnd2+13
                                    
                    if End2 < Posit < Strt2 :
                        AlignCount = AlignCount + 1
                        output.write(name+'\t'+str(tStrt2)+'\t'+str(tEnd2)+'\t'+str(Posit)+'\t'+count+'\n') #sense count
                
                  
            else :          # W orientation
                name = temp2[0]
                Orien2 = temp2[4]
                if Orien2 == '+' :  #also W orientation
                    tStrt2 = int(temp2[2])
                    tEnd2 = int(temp2[3])
                    Strt2 = tStrt2-16
                    End2 = tEnd2-13
                    if Strt2 < Posit < End2 :
                        AlignCount = AlignCount + 1
                        output.write(name+'\t'+str(tStrt2)+'\t'+str(tEnd2)+'\t'+str(Posit)+ '\t'+count+'\n') #sense count
                        AlignCountTotal = AlignCountTotal + float(count)
    file2.close()
print (AlignCountTotal)
output.close()