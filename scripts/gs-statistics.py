import statistics
import sys


def go(fname):
    inp = open(fname)
    l = []
    for line in inp:
        l.append(int(line.split('\t')[3]))
    inp.close()

    print('File: ', fname)
    print('Number: ', len(l))
    print('Min: ', min(l))
    print('Max: ', max(l))
    print('Mean: ', statistics.mean(l))
    print('Median: ', statistics.median(l))
    print('StdDev: ', statistics.stdev(l))
    print()


def main():
    go('/Users/alok/Tools/indel-detect/scripts/gs-svclassify.txt')
    go('/Users/alok/Tools/indel-detect/scripts/gs-dgv.txt')


if __name__ == '__main__':
    if sys.version_info[0] < 3:
        raise Exception('Use python3')

    main()
