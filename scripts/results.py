import subprocess

FLAG_SVCLASSIFY = True
FLAG_COV = False
FLAG_50 = False
FLAG_SUPPORT = False
VAL_LARGE_SUPPORT = 3
VAL_SMALL_SUPPORT = 3
ONLY_LARGE = False
ONLY_SMALL = False

# Reciprocal overlap
RO = 0.5
INP_COV = 30

# Chromosomes to verify
# verifyChr = ['18']
verifyChr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
             '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

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


def checkCoverage(x_list):
    x_chr = x_list[0]
    x_st = x_list[1]
    x_en = x_list[2]
    cmd = "samtools depth -aa -r " + str(x_chr) + ":" + str(x_st) + "-" + str(x_en) \
        + " /Users/alok/Data/30x/chr/" + str(x_chr) \
        + ".bam | awk '{ sum += $3; n++} END { if (n > 0) print sum/n; }'"
    cov = float(subprocess.check_output(
        cmd, shell=True).decode("utf-8") .replace('\n', ''))

    return cov <= INP_COV


def parse():
    if FLAG_SVCLASSIFY:
        ref = open('/Users/alok/Tools/indel-detect/scripts/gs-svclassify.txt')
    else:
        ref = open('/Users/alok/Tools/indel-detect/scripts/gs-dgv.txt')

    """ softsv """
    # inp = open('/Users/alok/Data/30x/results/softsv.txt')
    """ lumpy """
    # inp = open('/Users/alok/Data/30x/results/lumpy_dels.txt')
    """ my dels"""
    # inp = open('/Users/alok/tmp/Dec/18/18_filter_deletions.txt')
    # inp = open('/Users/alok/tmp/Dec/18/18.txt')
    inp = open('/Users/alok/tmp/Dec/18/dels.txt')
    ref_list = []
    num_pred = 0

    for ref_line in ref:
        ref_l = ref_line.split()
        if ref_l[0] not in verifyChr:
            continue
        ref_st = int(ref_l[1])
        ref_en = int(ref_l[2])
        ref_sz = ref_en - ref_st + 1
        if ONLY_LARGE:
            if ref_sz < 500:
                continue
        if ONLY_SMALL:
            if ref_sz >= 500:
                continue
        here_list = [modChr(ref_l[0]), ref_st, ref_en]
        ref_list.append(here_list)

    found_list = [[False, 0, 0]]*len(ref_list)

    for inp_line in inp:
        inp_l = inp_line.split()
        # if len(inp_l) < 5:
        #     continue
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
        if len(inp_l) >= 6:
            sup1 = int(inp_l[4])
            sup2 = int(inp_l[5])
            if FLAG_SUPPORT:
                if inp_en - inp_st < 500:
                    if sup2 < VAL_SMALL_SUPPORT:
                        continue
                else:
                    if sup1 + sup2 < VAL_LARGE_SUPPORT:
                        continue

        if FLAG_50:
            if inp_en - inp_st + 1 < 50:
                continue
        if FLAG_COV:
            if not checkCoverage(here_list):
                continue

        num_pred += 1
        # print(here_list, inp_en - inp_st + 1)

        for i in range(len(ref_list)):
            if checkOverlap(ref_list[i], here_list):
                found_list[i] = [True, sup1, sup2]

    ref.close()
    inp.close()

    num_ref = len(ref_list)
    num_inp = num_pred
    num_true = 0
    for i in range(len(found_list)):
        if found_list[i][0]:
            num_true += 1

    print('Ref Total: ' + str(num_ref))
    print('Inp Total: ' + str(num_inp))
    print('Found: ' + str(num_true))

    precision = float(num_true/num_inp)
    recall = float(num_true / num_ref)
    f_score = float(2*precision*recall/(precision+recall))

    print('Precision: %.3f' % precision)
    print('Recall: %.3f' % recall)
    print('F-score: %.3f' % f_score)

    large = 0
    small = 0
    large_co = 0
    small_co = 0
    for i in range(len(ref_list)):
        ref_st = int(ref_list[i][1])
        ref_en = int(ref_list[i][2])
        ref_len = ref_en - ref_st + 1
        if ref_len > 500:
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


def main():
    parse()

    # returns True
    # print(checkOverlap('1', '1', [100, 200], [130, 201]))
    # returns False
    # print(checkOverlap('1', '1', [100, 200], [180, 220]))


if __name__ == '__main__':
    main()
