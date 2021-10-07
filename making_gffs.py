fwdsam= open ('/Users/annemariestringer/Desktop/SRR8156054ret-trimmedrevcomp.sam','r')#use the sam mapped to the reverse genome here
fwdsam.readline()
fwdsam.readline()
fwdsam.readline()

genomeF = [0 for x in range (4558954)] #make a list of 4558954 entries, set them all to 0

column_4_fwd = [] 

for line in fwdsam:
    if line.split('\t')[1] == '16':
        column_4_fwd.append(4558953-int(line.split('\t')[3])) 
        

for entry in column_4_fwd:
    genomeF [entry] += 1

fwdsam.close()

revsam= open ('/Users/annemariestringer/Desktop/SRR8156054ret-trimmedfwd.sam','r')#use the sam mapped to the forward genome here
revsam.readline()
revsam.readline()
revsam.readline()

genomeR = [0 for x in range (4558953)] #make a list of 4558953 entries (len of genome +1), set them all to 0

column_4_rev = []

for line in revsam:
    if line.split('\t')[1] == '16': 
        column_4_rev.append(int(line.split('\t')[3])) 
                    
for entry in column_4_rev:
    genomeR [entry] += 1

revsam.close()

gff=open('/Users/annemariestringer/Desktop/ret_4_23_21.gff','w')

gff_label = ('ret') #define track label

for x in range (0,len(genomeF)):   #runs loop until the end of list (starts at beginning)
    if genomeF [x]>0: #dont want file full of 0s
        gff.write('NA\tAgilent\t')
        gff.write(gff_label)
        gff.write('\t')
        gff.write(str(x))
        gff.write('\t')
        gff.write(str(x))
        gff.write('\t')
        gff.write(str(genomeF[x]))
        gff.write('\t.\t.\t.\n')

for x in range (0,len(genomeR)):
    if genomeR [x]>0:
        gff.write('NA\tAgilent\t')
        gff.write(gff_label)
        gff.write('\t')
        gff.write(str(x))
        gff.write('\t')
        gff.write(str(x))
        gff.write('\t')
        gff.write(str(genomeR[x]*-1))
        gff.write('\t.\t.\t.\n')
gff.close()
                    




