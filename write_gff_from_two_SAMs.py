input_revcomp_SAM_file = input("Enter the path and filename for the reverse complement .SAM file\n")
input_fwd_SAM_file = input("Enter the path and filename for the forward .SAM file\n")
output_gff_file = input("Enter the path and filename for output .gff file\n")
output_gff_label = input("Enter the track label for the output .gff file\n")

revcompsam= open(input_revcomp_SAM_file,'r')
revcompsam.readline()
revcompsam.readline()
revcompsam.readline()

genomeF = [0 for x in range (4558954)] #make a list of 4558953 entries (len of genome +1), set them all to 0

column_4_rev = [] 

for line in revcompsam:
    if line.split('\t')[1] == '16':
        column_4_rev.append(4558953-int(line.split('\t')[3])) 
        

for entry in column_4_rev:
    genomeF [entry] += 1

revcompsam.close()

fwdsam= open(input_fwd_SAM_file, 'r')
fwdsam.readline()
fwdsam.readline()
fwdsam.readline()

genomeR = [0 for x in range (4558953)] #make a list of 4558953 entries (len of genome), set them all to 0

column_4_fwd = []

for line in fwdsam:
    if line.split('\t')[1] == '16': 
        column_4_fwd.append(int(line.split('\t')[3])) 
                    
for entry in column_4_fwd:
    genomeR [entry] += 1

fwdsam.close()


gff= open (output_gff_file, 'w')
gff_label = (output_gff_label) #define track label

for x in range (0,len(genomeF)): 
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
                    




