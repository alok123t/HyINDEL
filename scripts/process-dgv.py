inp = open('/Users/alok/Data/10x/GRCh37_hg19_variants_2014-10-16.txt')
out = open('/Users/alok/Tools/indel-detect/scripts/dgv.txt', 'w')

for line in inp:
    if line.find('NA12878') == -1:
        continue
    l = line.split()
    if l[0] == 'variantaccession':
        continue
    if l[5] != 'deletion':
        continue
    sz = int(l[3])- int(l[2]) + 1
    s = l[1] + '\t' + l[2] + '\t' + l[3] + '\t' + str(sz) + '\n'
    out.write(s)

inp.close()
out.close()