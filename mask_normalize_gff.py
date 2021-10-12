input_gff_file = input("Enter the path and filename for the .gff file to be modified\n")
output_gff_file = input("Enter the path and filename for the output .gff file\n")
mask_file = input("Enter the path and filename for the .txt file listing regions to mask\n")
output_gff_label = input("Enter the track label for the output .gff file\n")

with open(mask_file, 'r') as ncfasta:
    newgffstart = []
    newgffstop = []
    newgff = []
    a= ncfasta.readline()

    while len(a) > 2:
        if a[0] == '|':
            b = a.split('|')
            c = b[3].split(':')
            d= c[1].split('-')
            newgff.append(d[0] + '\t' + d[1] + '\n')
        if a[0] == '>':       
            b = a.split('|')
            c = b[2].split(':')
            d= c[1].split('-')
            newgff.append(d[0] + '\t' + d[1] + '\n')
        a= ncfasta.readline()  
    newgff.append('2621995' + '\t' + '2622357' + '\n')


with open(input_gff_file, 'r') as inputgff:
    rpmgff = []

    for b in inputgff:
        toggle = True
        for line in newgff:
            if int(line.split()[0]) <= int(b.split()[3]) <= int(line.split()[1]):
                toggle = False
                break
        if toggle == True:
            rpmgff.append(b)
gff_label = (output_gff_label) #define track label
with open(output_gff_file, 'w') as newrpmgff:
    covsum = 0
    for c in rpmgff:
        covsum = covsum + abs(int(c.split()[5]))
    covsum = covsum/1000000

    for c in rpmgff:
        newrpmgff.writelines(c.split()[0]+'\t'+c.split()[1]+'\t'+gff_label+'\t'+str(c.split()[3])+'\t'+str(c.split()[4])+'\t'+str(int(c.split()[5])/covsum)+'\t.\t.\t.\t' + '\n')
 
