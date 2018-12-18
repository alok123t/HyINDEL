inp = open(
    '/Users/alok/Tools/indel-detect/scripts/Personalis_1000_Genomes_deduplicated_deletions.bed')
out = open('/Users/alok/Tools/indel-detect/scripts/gs-svclassify.txt', 'w')

for l in inp:
    line = l.split('\t')
    if line[0] == 'Chr':
        continue
    l_chr = line[0]
    l_st = int(line[1])
    l_en = int(line[2].replace('\n', ''))
    sz = int(l_en) - int(l_st) + 1
    s = l_chr + '\t' + str(l_st) + '\t' + str(l_en) + '\t' + str(sz) + '\n'
    out.write(s)

inp.close()
out.close()
