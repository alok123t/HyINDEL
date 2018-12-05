inp = open('/Users/alok/tmp/29Nov/sim.bedpe')
out = open('/Users/alok/Tools/indel-detect/scripts/sim-events.txt', 'w')

for l in inp:
    line = l.split('\t')
    st = int(line[2])
    en = int(line[4])
    sz = en - st
    s = line[0] + '\t' + str(st) + '\t' + str(en) + '\t' + str(sz) + '\n'
    out.write(s)

inp.close()
out.close()
