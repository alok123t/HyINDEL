from random import choice, randint
from math import ceil
import sys

delsDistr = [(50, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 700), (700, 800), (800, 900),
             (900, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 5000), (5000, 6000), (6000, 7000), (7000, 8000), (8000, 9000), (9000, 10000)]
totalDels = 2000

insDistr = [(50, 100), (100, 150), (150, 200),
            (200, 250), (250, 300), (300, 350)]
totIns = 500


def main():
    if len(sys.argv) != 2:
        print('Usage: python range.py inpFile.bed')
        return

    inpName = sys.argv[1]
    inpFile = open(inpName, 'r')
    lines = inpFile.readlines()
    inpFile.close()

    delsCount = [0]*len(delsDistr)

    for line in lines:
        l = line.split('\t')
        if l[0] == 'Chr':
            continue
        st = int(l[1])
        en = int(l[2])
        sz = en - st + 1

        for i in range(len(delsDistr)):
            if delsDistr[i][0] <= sz and sz <= delsDistr[i][1]:
                delsCount[i] += 1
                break

    delsProb = [float(1.0 * x / sum(delsCount)) for x in delsCount]

    for i in range(len(delsDistr)):
        reps = int(ceil(totalDels * delsProb[i]))
        while reps:
            reps -= 1
            sz = randint(delsDistr[i][0], delsDistr[i][1])
            outStr = 'DEL ' + str(sz) + ' ' + str(sz) + ' 1'
            print(outStr)

    for i in range(len(insDistr)):
        reps = int(ceil(totIns / len(insDistr)))
        while reps:
            reps -= 1
            sz = randint(insDistr[i][0], insDistr[i][1])
            seq = ''.join(choice('ACGT') for _ in range(sz))
            outStr = 'INC ' + seq
            print(outStr)


if __name__ == '__main__':
    main()
