from random import choice
import sys

HETEROZYGOUS_VARS = 1000


def main():
    if len(sys.argv) != 2:
        print('Usage: python half.py pathToSimDir')
        return

    simDir = sys.argv[1]
    if simDir[-1] != '/':
        simDir += '/'

    inpBedpe = open(simDir + 'mod.bedpe', 'r')
    linesBedpe = inpBedpe.readlines()
    inpBedpe.close()

    out1 = open(simDir + 'sim1.txt', 'w')
    out2 = open(simDir + 'sim2.txt', 'w')

    co = 0
    for i in range(len(linesBedpe)):
        line = linesBedpe[i]
        l = line.split()
        if l[0] == 'LITERAL':
            sz = int(l[2])
            chr = l[3]
            pos = l[4]
            seq = ''.join(choice('ACGT') for _ in range(sz))
            out = 'INCL ' + chr + ' ' + \
                str(pos) + ' ' + str(pos) + ' ' + seq + '\n'
            out1.write(out)
            if co < HETEROZYGOUS_VARS:
                out2.write(out)
            co += 1
        elif 'DEL' in line:
            chr = l[0]
            st = int(l[2])
            en = int(l[4]) - 1
            out = 'DELL ' + chr + ' ' + str(st) + ' ' + str(en) + '\n'
            out1.write(out)
            if co < HETEROZYGOUS_VARS:
                out2.write(out)
            co += 1
        else:
            continue

    out1.close()
    out2.close()


if __name__ == '__main__':
    main()
