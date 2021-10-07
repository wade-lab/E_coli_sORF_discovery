with open('/Users/annemariestringer/Desktop/CP001509_ncRNA.fasta', 'r') as ncfasta:
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


with open('/Users/annemariestringer/Desktop/ret_4_23_21.gff', 'r') as retgff:
    rpmgff = []

    for b in retgff:
        toggle = True
        for line in newgff:
            if int(line.split()[0]) <= int(b.split()[3]) <= int(line.split()[1]):
                toggle = False
                break
        if toggle == True:
            rpmgff.append(b)

with open('/Users/annemariestringer/Desktop/4-23-21Ret-RPMgff.txt', 'w') as rpmretgff:
    covsum = 0
    for c in rpmgff:
        covsum = covsum + abs(int(c.split()[5]))
    covsum = covsum/1000000

    for c in rpmgff:
        rpmretgff.writelines(c.split()[0]+'\t'+c.split()[1]+'\t'+'Ret-RPM'+'\t'+str(c.split()[3])+'\t'+str(c.split()[4])+'\t'+str(int(c.split()[5])/covsum)+'\t.\t.\t.\t' + '\n')
 
"""   
    with open('/Users/carol/Desktop/Apiret project/Joe programs/Ret-RPMgff.txt', 'w') as rpmgff:

        for b in retgff:
            toggle = True
            for line in newgff:
                if int(line.split()[0]) <= int(b.split()[3]) <= int(line.split()[1]):
                    toggle = False
                    break
            if toggle == True:
                rpmgff.write(b)

                #rpmgff.write(b.split()[0]+'\t'+b.split()[1]+'\t'+'Ret-RPM'+'\t'+str(b.split()[3])+'\t'+str(b.split()[4])+'\t'+str(int(b.split()[5])/covsum)+'\t'+b.split()[6]+'\t'+b.split()[7]+'\t'+b.split()[8]+'\n')
"""        



"""
    b = newretgff.readline()   
    covsum = 0
    
    while len(b) > 2:
        covsum = covsum + abs(int(b.split()[5]))
        b = retgff.readline()
    covsum = covsum/1000000
#need to sum only lines from filtered gff
"""
