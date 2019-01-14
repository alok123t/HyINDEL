
FLAG_SVCLASSIFY = True
FLAG_50 = True
ONLY_LARGE = False
ONLY_SMALL = False
FLAG_MAPPABILITY = False
FILTER_MAPPABILITY = 0.75
SPLIT_LARGE = 500

# Reciprocal overlap
RO = 0.5

"""
This function normalizes chromosome name
chrchr1 -> 1
chr1 -> 1
1 -> 1
"""


def modChr(s):
    if s.find('chrchr') != -1:
        return s[6:]
    elif s.find('chr') != -1:
        return s[3:]
    else:
        return s


"""
This function returns
True: if there is reciprocal overlap and are on same chr
False: else
"""


def checkOverlap(list_a, list_b):
    if list_a[0] != list_b[0]:
        return False

    st1 = list_a[1]
    en1 = list_a[2]
    st2 = list_b[1]
    en2 = list_b[2]

    if (st1 <= st2 and st2 <= en1) or \
            (st2 <= st1 and st1 <= en2):
        st = max(st1, st2)
        en = min(en1, en2)
        overlap_len = float(en - st + 1)
        a_len = float(en1 - st1 + 1)
        b_len = float(en2 - st2 + 1)

        ro_a = overlap_len/a_len
        ro_b = overlap_len/b_len

        # print(ro_a, ro_b)

        if ro_a >= RO and ro_b >= RO:
            return True

    return False


def checkMap(st, en, sc):
    map_st = int(st / 100)
    map_en = int(en / 100)
    here = []
    for i in range(map_st-1, map_en):
        here.append(sc[i])
    map_sc = 1.0 * sum(here) / len(here)
    # print(st, ref_sz, map_sc, here)
    if map_sc < FILTER_MAPPABILITY:
        return True


def parse(fName, verifyChr):
    print(verifyChr, fName)
    if FLAG_SVCLASSIFY:
        ref = open(
            '/Users/alok/Tools/indel-detect/scripts/GS/gs-svclassify.txt')
    else:
        ref = open(
            '/Users/alok/Tools/indel-detect/scripts/GS/gs-dgv.txt')

    inp = open(fName)

    mapFile = open('/Users/alok/Data/mappability_hg19/chr18.dat.txt')
    lines_map = mapFile.readlines()
    sc = [float(x.strip()) for x in lines_map]

    ref_list = []
    num_pred = 0

    for ref_line in ref:
        ref_l = ref_line.split()
        if ref_l[0] not in verifyChr:
            continue
        ref_st = int(ref_l[1])
        ref_en = int(ref_l[2])
        ref_sz = ref_en - ref_st + 1
        if FLAG_MAPPABILITY:
            if checkMap(ref_st, ref_en, sc):
                continue
        if ONLY_LARGE:
            if ref_sz < SPLIT_LARGE:
                continue
        if ONLY_SMALL:
            if ref_sz >= SPLIT_LARGE:
                continue
        here_list = [modChr(ref_l[0]), ref_st, ref_en]
        ref_list.append(here_list)

    found_list = [[False, 0, 0]]*len(ref_list)

    for inp_line in inp:
        inp_l = inp_line.split()
        if inp_l[0] == 'Chromosome':
            continue
        inp_chr = modChr(inp_l[0])
        if inp_chr not in verifyChr:
            continue
        inp_st = int(inp_l[1])
        inp_en = int(inp_l[2])
        here_list = [inp_chr, inp_st, inp_en]
        sup1 = 0
        sup2 = 0
        if FLAG_MAPPABILITY:
            if checkMap(inp_st, inp_en, sc):
                continue
        if ONLY_LARGE:
            if inp_en - inp_st + 1 < SPLIT_LARGE:
                continue
        if ONLY_SMALL:
            if inp_en - inp_st + 1 >= SPLIT_LARGE:
                continue

        if FLAG_50:
            if inp_en - inp_st + 1 < 50:
                continue

        num_pred += 1
        # print(here_list, inp_en - inp_st + 1)

        for i in range(len(ref_list)):
            if checkOverlap(ref_list[i], here_list):
                found_list[i] = [True, sup1, sup2]

    ref.close()
    inp.close()
    mapFile.close()

    num_ref = len(ref_list)
    num_inp = num_pred
    num_true = 0
    for i in range(len(found_list)):
        if found_list[i][0]:
            num_true += 1

    print('Ref Total: ' + str(num_ref))
    print('Inp Total: ' + str(num_inp))
    print('Found: ' + str(num_true))

    precision = float(1.0*num_true/num_inp)
    recall = float(1.0*num_true / num_ref)
    print('Precision: %.3f' % precision)
    print('Recall: %.3f' % recall)
    f_score = float(2.0*precision*recall/(precision+recall))
    print('F-score: %.3f' % f_score)

    large = 0
    small = 0
    large_co = 0
    small_co = 0
    for i in range(len(ref_list)):
        ref_st = int(ref_list[i][1])
        ref_en = int(ref_list[i][2])
        ref_len = ref_en - ref_st + 1
        if ref_len > SPLIT_LARGE:
            large += 1
            if found_list[i][0]:
                large_co += 1
            # else:
            #     print('Large:', ref_st, ref_en,
            #           found_list[i][1], found_list[i][2])
            # print('chr' + verifyChr[0] + ":" + str(ref_st) + "-" +
            #       str(ref_en), ref_len, found_list[i][0])
        else:
            small += 1
            if found_list[i][0]:
                small_co += 1
            # else:
            #     print('Small:', ref_st, ref_en,
            #           found_list[i][1], found_list[i][2])
            # print('chr' + verifyChr[0] + ":" + str(ref_st) + "-" +
            #       str(ref_en), ref_len, found_list[i][0])
    print(large, small)
    print(large_co, small_co)
    print('---------------')


def main():
    f_softsv = '/Users/alok/Data/30x/results/softsv.txt'
    # f_lumpy = '/Users/alok/tmp/2019/Jan/4/dels-lumpy.txt'
    f_split = '/Users/alok/tmp/2018/Dec/20/all_split.bed'
    f_merge = '/Users/alok/tmp/2018/Dec/20/merge_merge.bed'
    f_sc = '/Users/alok/Data/30x/results/all.bed'
    f_tiddit = '/Users/alok/Data/30x/results/tiddit/noins.txt'
    f_small_mq = '/Users/alok/tmp/2019/Jan/14/1/tmp1_dels_small.txt'
    f_large_mq = '/Users/alok/tmp/2019/Jan/14/18/tmp1_dels_large.txt'
    f_all_mq = '/Users/alok/tmp/2019/Jan/14/19/deletions.bed'
    f_mq_18 = '/Users/alok/tmp/2019/Jan/12/18_mq/deletions.bed'
    f_manta = '/Users/alok/Data/30x/results/manta.txt'

    allChr = [['19']]
    # allChr = [['1'], ['2'], ['3'], ['4'], ['5'], ['6'], ['7'], ['8'], ['9'], ['10'], ['11'], [
    #     '12'], ['13'], ['14'], ['15'], ['16'], ['17'], ['18'], ['19'], ['20'], ['21'], ['22'], ['X']]

    # for curChr in allChr:
    #     # parse(f_softsv, curChr)
    #     parse(f_lumpy, curChr)
    #     parse(f_our, curChr)

    wholeGenom = [x for sublist in allChr for x in sublist]
    # parse(f_softsv, wholeGenom)
    # parse(f_lumpy, wholeGenom)
    # parse(f_sc, wholeGenom)
    # parse(f_split, wholeGenom)
    parse(f_merge, wholeGenom)
    parse(f_all_mq, wholeGenom)
    # parse(f_large_mq, wholeGenom)
    # parse(f_small_mq, wholeGenom)
    # parse(f_mq_18, wholeGenom)
    # parse(f_tiddit, wholeGenom)
    # parse(f_manta, wholeGenom)

    # returns True
    # print(checkOverlap('1', '1', [100, 200], [130, 201]))
    # returns False
    # print(checkOverlap('1', '1', [100, 200], [180, 220]))


if __name__ == '__main__':
    main()
