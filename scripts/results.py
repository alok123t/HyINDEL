
FLAG_SVCLASSIFY = 0
FLAG_50 = True
ONLY_LARGE = False
ONLY_SMALL = False
FLAG_MAPPABILITY = False
FILTER_MAPPABILITY = 0.75
SPLIT_LARGE = 500
DIST_INS = 5

PRINT_MISSING = False

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


def checkOverlapDels(list_a, list_b):
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

        if a_len == 0 or b_len == 0:
            return False

        ro_a = overlap_len/a_len
        ro_b = overlap_len/b_len

        if ro_a >= RO and ro_b >= RO:
            return True

    return False


def checkOverlapIns(list_a, list_b):
    if list_a[0] != list_b[0]:
        return False
    return abs(list_a[1] - list_b[1]) <= DIST_INS


def checkMap(st, en, sc):
    map_st = int(st / 100)
    map_en = int(en / 100)
    here = []
    for i in range(map_st-1, map_en):
        here.append(sc[i])
    map_sc = 1.0 * sum(here) / len(here)
    if map_sc < FILTER_MAPPABILITY:
        return True


def readMy(fName, delsFlag):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if delsFlag:
            chr = l[0]
            st = int(l[1])
            en = int(l[2])
            here = [modChr(chr), st, en]
            ret.append(here)
        else:
            chr = l[0]
            st = int(l[1])
            here = [modChr(chr), st]
            ret.append(here)

    return ret


def readLumpy(fName):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if l[4] != '<DEL>':
            continue
        chr = l[0]
        st = int(l[1])
        info = l[7].split(';')
        en = int(info[3].split('=')[1])
        here = [modChr(chr), st, en]
        ret.append(here)

    return ret


def readTiddit(fName):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if l[4] != '<DEL>':
            continue
        if l[6] != 'PASS':
            continue
        chr = l[0]
        st = int(l[1])
        info = l[7].split(';')
        en = int(info[3].split('=')[1])
        here = [modChr(chr), st, en]
        ret.append(here)

    return ret


def readGrom(fName, delsFlag):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if delsFlag:
            if l[4] != '<DEL>':
                continue
            chr = l[0]
            st = int(l[1])
            info = l[7].split(';')
            en = int(info[1].split('=')[1])
            here = [modChr(chr), st, en]
            ret.append(here)
        else:
            if l[4] != '<INS>':
                continue
            chr = l[0]
            st = int(l[1])
            here = [modChr(chr), st]
            ret.append(here)

    return ret


def readSoftsv(fNameSmall, fNameLarge):
    inp = open(fNameSmall)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'Chromosome':
            continue
        chr = l[0]
        st = int(l[1])
        en = int(l[2])

        here = [modChr(chr), st, en]
        ret.append(here)

    inp = open(fNameLarge)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'Chromosome':
            continue
        chr = l[0]
        st = int(l[1])
        en = int(l[2])
        here = [modChr(chr), st, en]
        ret.append(here)

    return ret


def readSvclassify(fName):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        chr = l[0]
        if chr == 'Chr':
            continue
        st = int(l[1])
        en = int(l[2])

        here = [modChr(chr), st, en]
        ret.append(here)

    return ret


def readSim(fName, delsFlag):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'LITERAL':
            continue
        if delsFlag:
            if 'DEL' in l[6]:
                chr = l[0]
                st = int(l[2])
                en = int(l[4])
                here = [modChr(chr), st, en]
                ret.append(here)
        else:
            if 'INC' in l[6]:
                chr = l[0]
                st = int(l[2])
                here = [chr, st]
                ret.append(here)

    return ret


def compare_dels(tool, ref, toolName):
    found = [False] * len(ref)

    for i in range(len(tool)):
        for j in range(len(ref)):
            if checkOverlapDels(tool[i], ref[j]):
                found[j] = True

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = sum(found)

    precision = float(1.0 * num_true / num_tool)
    recall = float(1.0 * num_true / num_ref)

    print(toolName)
    print('Precision: %.3f' % precision)
    print('Recall: %.3f' % recall)
    f_score = float(2.0*precision*recall/(precision+recall))
    print('F-score: %.3f' % f_score)
    print('--------------------')


def compare_ins(tool, ref, toolName):
    found = [False] * len(ref)

    for i in range(len(tool)):
        for j in range(len(ref)):
            if checkOverlapIns(tool[i], ref[j]):
                found[j] = True

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = sum(found)

    precision = float(1.0 * num_true / num_tool)
    recall = float(1.0 * num_true / num_ref)

    print(toolName)
    print('Precision: %.3f' % precision)
    print('Recall: %.3f' % recall)
    f_score = float(2.0*precision*recall/(precision+recall))
    print('F-score: %.3f' % f_score)
    print('--------------------')


def parse(fName, verifyChr):
    print(verifyChr, fName)
    if FLAG_SVCLASSIFY == 0:
        ref = open(
            '/Users/alok/Tools/indel-detect/scripts/GS/gs-svclassify.txt')
    elif FLAG_SVCLASSIFY == 1:
        ref = open(
            '/Users/alok/tmp/today/mills.bed')
    elif FLAG_SVCLASSIFY == 2:
        ref = open(
            '/Users/alok/tmp/today/long.bed')

    inp = open(fName)

    mapFile = open('/Users/alok/Data/mappability_hg19/chr18.dat.txt')
    lines_map = mapFile.readlines()
    sc = [float(x.strip()) for x in lines_map]

    ref_list = []
    num_pred = 0

    for ref_line in ref:
        ref_l = ref_line.split()
        ref_chr = modChr(ref_l[0])
        if ref_chr not in verifyChr:
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
        found = False
        for i in range(len(ref_list)):
            if checkOverlap(ref_list[i], here_list):
                found_list[i] = [True, sup1, sup2]
                found = True
        # if not found:
        #     print(here_list[0] + ':' + str(here_list[1])+'-'+str(here_list[2]))

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
            else:
                if PRINT_MISSING:
                    print('chr' + ref_list[i][0] + ":" + str(ref_st-450) + "-" +
                          str(ref_en+450))
                #     print('Large:', ref_st, ref_en,
                #           found_list[i][1], found_list[i][2])

        else:
            small += 1
            if found_list[i][0]:
                small_co += 1
            else:
                if PRINT_MISSING:
                    print(ref_list[i][0] + ":" + str(ref_st) + "-" +
                          str(ref_en))
                # print('Small:', ref_st, ref_en,
                #       found_list[i][1], found_list[i][2])
    print(large, small)
    print(large_co, small_co)
    print('---------------')


def main():
    real_lumpy = readLumpy(
        '/Users/alok/tmp/2019/Feb/19/Results/Real/lumpy_real.vcf')
    real_tiddit = readTiddit(
        '/Users/alok/tmp/2019/Feb/19/Results/Real/tiddit_real.vcf')
    real_svclassify_ref = readSvclassify(
        '/Users/alok/Tools/indel-detect/scripts/GS/Personalis_1000_Genomes_deduplicated_deletions.bed')

    sim_dels_lumpy = readLumpy('/Users/alok/tmp/today/out_rg.vcf')
    sim_dels_tiddit = readTiddit('/Users/alok/tmp/today/out.vcf')
    sim_dels_softsv = readSoftsv(
        '/Users/alok/tmp/today/deletions_small.txt', '/Users/alok/tmp/today/deletions.txt')
    sim_dels_grom = readGrom('/Users/alok/tmp/today/out.SV.vcf', True)
    sim_ins_grom = readGrom('/Users/alok/tmp/today/out.SV.vcf', False)
    sim_dels_my = readMy(
        '/Users/alok/tmp/2019/Feb/19/Results/Simulations/30x/deletions.bed', True)
    sim_ins_my = readMy('/Users/alok/tmp/2019/Feb/22/30x/tmp1_ins.txt', False)
    sim_dels_ref = readSim(
        '/Users/alok/tmp/2019/Feb/14/Simulations/sim_ref.bedpe', True)
    sim_ins_ref = readSim(
        '/Users/alok/tmp/2019/Feb/14/Simulations/sim_ref.bedpe', False)

    compare_dels(real_lumpy, real_svclassify_ref, 'REAL-DELS-LUMPY')
    compare_dels(real_tiddit, real_svclassify_ref, 'REAL-DELS-TIDDIT')

    compare_dels(sim_dels_lumpy, sim_dels_ref, 'SIM-DELS-LUMPY')
    compare_dels(sim_dels_tiddit, sim_dels_ref, 'SIM-DELS-TIDDIT')
    compare_dels(sim_dels_softsv, sim_dels_ref, 'SIM-DELS-SOFTSV')
    compare_dels(sim_dels_grom, sim_dels_ref, 'SIM-DELS-GROM')
    compare_dels(sim_dels_my, sim_dels_ref, 'SIM-DELS-MY')

    compare_ins(sim_ins_grom, sim_ins_ref, 'SIM-INS-GROM')
    compare_ins(sim_ins_my, sim_ins_ref, 'SIM-INS-MY')

    f_softsv = '/Users/alok/Data/30x/results/softsv.txt'
    f_lumpy = '/Users/alok/tmp/today/lumpy_dels.bed'
    f_tiddit = '/Users/alok/tmp/today/tiddit_dels.bed'
    f_manta = '/Users/alok/Data/30x/results/manta.txt'
    f_split = '/Users/alok/tmp/2018/Dec/20/all_split.bed'
    f_sc = '/Users/alok/Data/30x/results/all.bed'
    f_split_mq = '/Users/alok/tmp/2019/Jan/17/full/tmp2_split.bed'
    f_1 = '/Users/alok/tmp/2018/Dec/20/merge_merge.bed'
    f_2 = '/Users/alok/tmp/2019/Jan/17/full/deletions.bed'
    f_3 = '/Users/alok/tmp/2019/Jan/22/full/deletions.bed'
    f_3_1 = '/Users/alok/tmp/2019/Feb/8/windows_10/deletions.bed'
    f_3_2 = '/Users/alok/tmp/today/Results/Real/my/deletions.bed'
    f_3_3 = '/Users/alok/tmp/today/Results/Real/my_filter/deletions.bed'

    # allChr = [['18']]
    allChr = [['1'], ['2'], ['3'], ['4'], ['5'], ['6'], ['7'], ['8'], ['9'], ['10'], ['11'], [
        '12'], ['13'], ['14'], ['15'], ['16'], ['17'], ['18'], ['19'], ['20'], ['21'], ['22'], ['X']]

    # for curChr in allChr:
    #     parse(f_2, curChr)

    wholeGenom = [x for sublist in allChr for x in sublist]
    # parse(f_lumpy, wholeGenom)
    # parse(f_tiddit, wholeGenom)
    # parse(f_softsv, wholeGenom)
    # parse(f_manta, wholeGenom)
    # parse(f_sc, wholeGenom)
    # parse(f_split, wholeGenom)
    # parse(f_split_mq, wholeGenom)
    # parse(f_1, wholeGenom)
    # parse(f_2, wholeGenom)
    # parse(f_3, wholeGenom)
    # parse(f_3_1, wholeGenom)
    # parse(f_3_2, wholeGenom)
    # parse(f_3_3, wholeGenom)


if __name__ == '__main__':
    main()
