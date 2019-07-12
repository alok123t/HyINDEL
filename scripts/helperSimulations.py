import random
import sys

# Small, large size ranges
DISTR = [[50, 500], [500, 10000]]
# Number of variants for each size range
N_DISTR = 375


def getRandomString(l):
    return ''.join(random.choice('ATGC') for _ in range(l))


def genPositions():
    # For each size range
    for d in DISTR:
        # Insertions
        rndVals = random.sample(range(d[0], d[1]), N_DISTR)
        for val in rndVals:
            out = 'INC ' + getRandomString(val)
            print(out)
        # Deletions
        rndVals = random.sample(range(d[0], d[1]), N_DISTR)
        for val in rndVals:
            out = 'DEL ' + str(val)
            print(out)


def genVariants(folderPath):
    dels_large = []
    dels_small = []
    ins_large = []
    ins_small = []

    # Read input bedpe file
    with open(folderPath + 'tmp.bedpe', 'r') as f:
        for line in f:
            l = line.split('\t')
            if 'DEL' in l[6]:
                chr = l[0]
                st = int(l[2])
                en = int(l[4])
                if en - st + 1 < 500:
                    dels_small.append([chr, st, en])
                else:
                    dels_large.append([chr, st, en])
            elif 'LITERAL' in l[0]:
                chr = l[3]
                pos = int(l[4])
                sz = int(l[1])
                seq = getRandomString(sz)
                if sz < 500:
                    ins_small.append([chr, pos, seq])
                else:
                    ins_large.append([chr, pos, seq])

    # Write variants
    with open(folderPath + '1.txt', 'w') as f1, open(folderPath + '2.txt', 'w') as f2:
        # Dels large
        for i in range(len(dels_large)):
            x = dels_large[i]
            out = 'DELL ' + x[0] + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'
            if i < len(dels_large) / 2:
                f1.write(out)
            f2.write(out)
        # Dels small
        for i in range(len(dels_small)):
            x = dels_small[i]
            out = 'DELL ' + x[0] + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'
            if i < len(dels_small) / 2:
                f1.write(out)
            f2.write(out)
        # Ins large
        for i in range(len(ins_large)):
            x = ins_large[i]
            out = 'INCL ' + x[0] + ' ' + \
                str(x[1]) + ' ' + str(x[1]) + ' ' + x[2] + '\n'
            if i < len(ins_large) / 2:
                f1.write(out)
            f2.write(out)
        # Ins small
        for i in range(len(ins_small)):
            x = ins_small[i]
            out = 'INCL ' + x[0] + ' ' + \
                str(x[1]) + ' ' + str(x[1]) + ' ' + x[2] + '\n'
            if i < len(ins_small) / 2:
                f1.write(out)
            f2.write(out)


def main():

    if sys.argv[1] == 'genPositions':
        genPositions()
    elif sys.argv[1] == 'genVariants':
        genVariants(sys.argv[2])


if __name__ == '__main__':
    main()
