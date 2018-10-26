import os

# f1 = open('/Users/alok/tmp/my_del')
f1 = open('/Users/alok/tmp/1')
f2 = open('/Users/alok/Tools/indel-detect/scripts/dgv.bed')
done = ['chr1']
# done = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

def mod(x):
    if x.find('chrchr') != -1:
        return x[6:]
    elif x.find('chr') != -1:
        return x[3:]
    else:
        return x


def overlap(c_a, c_b, a, b):
    if mod(c_a) != mod(c_b):
        return False
    st1 = a[0]
    en1 = a[1]
    st2 = b[0]
    en2 = b[1]

    if (st2 >= st1 and st2 <= en1) or (st1 >= st2 and st1 <= en2):
        x = max(st1, st2)
        y = min(en1, en2)
        l = (y - x + 1)*1.0
        la = en1 - st1 + 1
        lb = en2 - st2 + 1
        # print(l, la, lb, l/la, l/lb)
        if l / la >= 0.5 and l / lb >= 0.5:
            return True    
    return False

# print(overlap("chr1", "chr1", [100, 200], [130, 201]))
# print(overlap("chr1", "chr1", [100, 200], [180, 220]))

chrz = []
res = []
for line in f2:
    l = line.split()
    if l[0] == 'Chromosome':
        continue
    if l[0] not in done:
        continue
    # if int(l[3]) <= 10000:
    if int(l[3]) > 10000 or int(l[3]) < 200:
        continue
    chrz.append(l[0])
    res.append([int(l[1]), int(l[2])])
print('Total : ' + str(len(res)))
found = [0]*len(res)
pred = 0
for line in f1:
    l = line.split()
    if len(l) < 7:
        continue
    if l[0] == 'Chromosome':
        continue
    if int(l[3]) > 10000 or int(l[3]) < 200:
        continue
    here = [int(l[1]), int(l[2])]
    # if int(l[4]) < 2 or int(l[5]) < 2:
    #     continue
    pred += 1
    for i in range(len(res)):
        if overlap(chrz[i], l[0], res[i], here):
            found[i] = 1
            # print(l[4] + ' ' + l[5])
            break
        #no break?
print('Found : ' + str(found.count(1)))
print('Total predictions : ' + str(pred))
for i in range(len(res)):
    x = res[i][0]
    y = res[i][1]
    cmd = "samtools depth -aa -r 1:" + str(x) + "-" + str(y) + " /Users/alok/Data/30x/chr1.bam     | awk '{ sum += $3; n++} END { if (n > 0) print sum/n; }'"
    print (found[i], end=' ', flush=True)
    os.system(cmd)
# for i in range(len(found)):
#     if found[i] > 1:
#         print(res[i])
#         break
