from collections import Counter


def check_overlap(fName):
    inp = open(fName)

    ar = []
    for line in inp:
        l = line.replace('\n', '').split('\t')
        ar.append(l)

    inp.close()

    f = []
    mp = {}

    for i in range(len(ar)):
        chr_1 = ar[i][0]
        st_1 = int(ar[i][1])
        en_1 = int(ar[i][2])
        len_1 = int(ar[i][3])

        for j in range(i + 1, len(ar)):
            chr_2 = ar[j][0]
            st_2 = int(ar[j][1])
            en_2 = int(ar[j][2])
            len_2 = int(ar[j][3])

            if chr_1 != chr_2:
                continue

            if st_1 <= st_2 and st_2 <= en_1:
                f.append([i, j])
                if i not in mp:
                    mp[i] = 0
                mp[i] += 1
                # print(i, j, len_1, len_2)

    print(max(mp.values()))
    print(Counter(mp.values()))


check_overlap('/Users/alok/Tools/indel-detect/scripts/gs-dgv.txt')
check_overlap('/Users/alok/Tools/indel-detect/scripts/svclassify.txt')
