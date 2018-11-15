inp = open('/Users/alok/Data/GS/NA12878_svs.vcf')
out = open('/Users/alok/Tools/indel-detect/scripts/gs-svclassify.txt', 'w')

for l in inp:
    if l[0] == '#':
        continue
    line = l.split('\t')
    if line[6] != "PASS":
        continue
    if line[4] != "<DEL>":
        continue
    st = line[1]
    en = line[7].split(';')[0].split('=')[1]
    sz = int(en) - int(st) + 1
    s = line[0] + '\t' + st + '\t' + en + '\t' + str(sz) + '\n'
    out.write(s)

inp.close()
out.close()
