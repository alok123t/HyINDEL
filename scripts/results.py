# pylint: disable=unused-variable

from statistics import median
import numpy as np
import matplotlib.pyplot as plt

# Remove events < 50 bp in tools
FLAG_50 = True
# Plot breakpoint error plot
BREAKPOINT_PLOT = False
# Plot support
SUPPORT_PLOT = True
# Print metrics seperately for large and small variants
PRINT_SPLIT = False
SPLIT_LARGE = 500
# Maximum error in breakpoint for insertions
DIST_INS = 10
# Reciprocal overlap
RO = 0.5
GT_METRIC = False

"""
This function normalizes chromosome name
chrchr1 -> 1
chr1 -> 1
1 -> 1
"""


class Deletion:
    def __init__(self, t_chr, t_st, t_en, t_iD=-1, t_pe=-1, t_sc=-1, t_sr=-1, t_isHomo="N/A", t_isSmall="N/A"):
        self.chr = modChr(t_chr)
        self.st = t_st
        self.en = t_en
        self.iD = t_iD
        self.pe = t_pe
        self.sc = t_sc
        self.sr = t_sr
        self.isHomo = t_isHomo
        self.isSmall = t_isSmall

    def __str__(self):
        return self.chr+'\t' + str(self.st)+'\t' + str(self.en)


class Insertion:
    def __init__(self, t_chr, t_pos, t_seq='', t_seqLen=-1, t_refPos=-1, t_iD=-1, t_sc=-1, t_sup=-1, t_isHomo="N/A", t_isSmall="N/A"):
        self.chr = modChr(t_chr)
        self.pos = t_pos  # Position on reference
        self.seq = t_seq
        self.seqLen = t_seqLen
        self.refPos = t_refPos  # Position on sample
        self.iD = t_iD
        self.sc = t_sc
        self.sup = t_sup
        self.isHomo = t_isHomo
        self.isSmall = t_isSmall

    def __str__(self):
        if len(self.seq) == 0:
            return self.chr + '\t' + str(self.pos)
        else:
            return self.chr + '\t' + str(self.pos) + '\t' + str(len(self.seq))


def modChr(s):
    if s.find('chrchr') != -1:
        return s[6:]
    elif s.find('chr') != -1:
        return s[3:]
    else:
        return s


def checkDel(a, b):
    if a.chr != b.chr:
        return False

    a_len = a.en - a.st + 1
    b_len = b.en - b.st + 1

    o_st = max(a.st, b.st)
    o_en = min(a.en, b.en)
    o_len = o_en - o_st + 1

    if a_len == 0 or b_len == 0:
        return False

    ro_a = o_len / a_len
    ro_b = o_len / b_len

    if ro_a >= RO and ro_b >= RO:
        return True
    else:
        return False


def checkIns(a, b):
    if a.chr != b.chr:
        return False

    f_pos = False
    if abs(a.pos - b.pos) <= DIST_INS:
        f_pos = True

    return f_pos


def readHyINDEL(fName, delsFlag):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split()
            if delsFlag:
                if l[4] != '<DEL>':
                    continue
                info_col = l[7].split(';')
                if l[6] != 'PASS':
                    continue
                format_col = l[9].split(':')
                gt = format_col[0]
                isHomo = "N/A"
                if gt == '1/1':
                    isHomo = True
                elif gt == '0/1':
                    isHomo = False
                sup_pe = int(format_col[2])
                sup_sr = int(format_col[3])
                sup_sc = int(format_col[4])

                ret.append(Deletion(t_chr=l[0], t_st=int(
                    l[1]), t_en=int(info_col[2].split('=')[1]), t_pe=sup_pe, t_sc=sup_sc, t_sr=sup_sr, t_isHomo=isHomo))
            else:
                if l[4] != '<INS>':
                    continue
                if 'IMPRECISE' in l[7]:
                    ret.append(Insertion(t_chr=l[0], t_pos=int(l[1])))
                else:
                    info_col = l[7].split(';')
                    seq = info_col[3].split('=')[1]
                    ret.append(Insertion(t_chr=l[0], t_pos=int(
                        l[1]), t_seq=seq, t_seqLen=len(seq)))

    return ret


def readLumpy(fName):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split()
            if l[4] != '<DEL>':
                continue
            info_col = l[7].split(';')
            format_col = l[9].split(':')
            sup_pe = int(format_col[2])
            sup_sr = int(format_col[3])
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                info_col[3].split('=')[1]), t_pe=sup_pe, t_sr=sup_sr))

    return ret


def readTiddit(fName):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split()
            if l[4] != '<DEL>':
                continue
            # if l[6] != 'PASS':
            #     continue
            info_col = l[7].split(';')
            format_col = l[9].split(':')
            sup_pe = int(format_col[2])
            sup_sr = int(format_col[3])
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                info_col[3].split('=')[1]), t_pe=sup_pe, t_sr=sup_sr))

    return ret


def readSoftsv(fNameSmall, fNameLarge):
    ret = []
    with open(fNameLarge, 'r') as fLarge, open(fNameSmall, 'r') as fSmall:
        for line in fSmall:
            l = line.split()
            if l[0] == 'Chromosome':
                continue
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                l[2]), t_pe=int(l[4]), t_sc=int(l[5])))

        for line in fLarge:
            l = line.split()
            if l[0] == 'Chromosome':
                continue
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                l[2]), t_pe=int(l[4]), t_sc=int(l[5])))
    return ret


def readPamir(fName):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if l[4] != '<INS>':
            continue
        if l[6] != 'PASS':
            continue
        info_col = l[7].split(';')
        length = int(info_col[1].split('=')[1])
        if length < 50:
            continue
        sup = int(info_col[4].split('=')[1])
        seq = info_col[5].split('=')[1]
        ret.append(Insertion(t_chr=l[0], t_pos=int(
            l[1]), t_seq=seq, t_sup=sup))

    return ret


def readPopins(fName):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split('\t')
            if l[7] == 'NOANCHOR':
                continue
            f = False
            for x in ret:
                if x.chr == l[0] and abs(x.pos - int(l[1])) < 5:
                    f = True
                    break
            if f:
                continue
            ret.append(Insertion(t_chr=l[0], t_pos=int(l[1])))

    return ret


def readSim(fName2, fName1, delsFlag):
    ret = []
    with open(fName2, 'r') as f2, open(fName1, 'r') as f1:
        d = {}
        for line in f1:
            d[line] = 1
        for line in f2:
            isHomo = False
            if line in d:
                isHomo = True
            l = line.split()
            if delsFlag:
                if l[0] == 'DELL':
                    en = int(l[3])
                    st = int(l[2])
                    sz = en - st + 1
                    isSmall = False
                    if sz <= 500:
                        isSmall = True
                    ret.append(
                        Deletion(t_chr=l[1], t_st=st, t_en=en, t_isHomo=isHomo, t_isSmall=isSmall))
            else:
                if l[0] == 'INCL':
                    seq = l[4]
                    seqLen = len(l[4])
                    isSmall = False
                    if seqLen <= 500:
                        isSmall = True
                    ret.append(Insertion(t_chr=l[1], t_pos=int(
                        l[2]), t_seq=seq, t_seqLen=seqLen, t_isHomo=isHomo, t_isSmall=isSmall))
    return ret


def readSvclassify(fName, delsFlag):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'Chr':
            continue
        if delsFlag:
            st = int(l[1])
            en = int(l[2])
            sz = en - st + 1
            isSmall = False
            if sz <= 500:
                isSmall = True
            ret.append(Deletion(t_chr=l[0], t_st=st,
                                t_en=en, t_isSmall=isSmall))
        else:
            ret.append(Insertion(t_chr=l[0], t_pos=int(l[1])))

    return ret


def readDgvNovelIns(fName):
    ret = []
    ret_mob = []
    ret_ins = []
    with open(fName, 'r') as f:
        for line in f:
            l = line.split('\t')
            if l[0] == 'variantaccession':
                continue
            if 'NA12878' not in line:
                continue
            if l[6] != '1000_Genomes_Consortium_Pilot_Project':
                continue
            if l[5] == 'novel sequence insertion':
                ret.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
            if l[5] == 'mobile element insertion':
                ret_mob.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
            if l[5] == 'insertion':
                ret_ins.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
    return ret, ret_mob, ret_ins


def readMergedBenchmark(fName, delsFlag):
    ret = []
    f = open(fName)

    for line in f:
        l = line.split('\t')
        infoCol = l[7]
        if l[6] != 'PASS':
            continue
        if delsFlag:
            if l[2] == 'DEL':
                st = int(l[1])
                en = st + int(infoCol.split(';')[1].split('=')[1]) - 1
                sz = en - st + 1
                isSmall = False
                if sz <= 500:
                    isSmall = True
                ret.append(
                    Deletion(t_chr=l[0], t_st=st, t_en=en, t_isSmall=isSmall))
        else:
            if l[2] == 'INS':
                le = int(infoCol.split(';')[1].split('=')[1])
                ret.append(Insertion(t_chr=l[0], t_pos=int(l[1]), t_seqLen=le))

    f.close()
    return ret


"""
returns [[breakpoint error], [found], [support], [length error]]
[breakpoint error]: only for detected variants
[found]: True/False if present in benchmark
[support]: sc/sr support only for detected variants (if present)
[length error]: only for precise insertions detected variants
"""


def compare(toolUnfiltered, ref, toolName, delsFlag, checkHomoAndSmall=False, checkSmall=False):
    print(toolName)
    tool = []
    if delsFlag:
        if FLAG_50:
            tool = list(
                filter(lambda x: (x.en - x.st + 1) >= 50, toolUnfiltered))
        else:
            tool = toolUnfiltered
    else:
        tool = toolUnfiltered

    found = [False] * len(ref)

    dis = []
    sup = []
    le = []
    gt_true = 0
    gt_false = 0
    for i in range(len(tool)):
        for j in range(len(ref)):
            if delsFlag:
                if checkDel(tool[i], ref[j]):
                    found[j] = True
                    dis.append(abs(tool[i].en - ref[j].en) +
                               abs(tool[i].st - ref[j].st))
                    exact_sup = tool[i].sc + tool[i].sr
                    if GT_METRIC:
                        if tool[i].isHomo == ref[j].isHomo:
                            gt_true += 1
                        else:
                            gt_false += 1
                    if exact_sup > 0:
                        sup.append(exact_sup)
            else:
                if checkIns(tool[i], ref[j]):
                    found[j] = True
                    dis.append(abs(tool[i].pos - ref[j].pos))
                    if tool[i].seqLen != -1:
                        le.append(abs(tool[i].seqLen - ref[j].seqLen))

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = sum(found)

    precision = float(100 * num_true / num_tool)
    recall = float(100 * num_true / num_ref)
    fScore = float(2 * precision*recall/(precision+recall))
    print('P: %.3f R: %.3f F: %.3f' % (precision, recall, fScore))

    # For real data
    if checkSmall:
        ref_small = 0
        ref_large = 0
        small = 0
        large = 0
        for i in range(len(ref)):
            if ref[i].isSmall == True:
                ref_small += 1
                if found[i]:
                    small += 1
            elif ref[i].isSmall == False:
                ref_large += 1
                if found[i]:
                    large += 1
        print('Small: %.3f Large: %.3f' %
              (small / ref_small * 100, large / ref_large * 100))

    # For simulated data
    if checkHomoAndSmall:
        homo_small = 0
        homo_large = 0
        ref_homo_small = 0
        ref_homo_large = 0
        hetero_small = 0
        hetero_large = 0
        ref_hetero_small = 0
        ref_hetero_large = 0
        for i in range(len(ref)):
            if ref[i].isHomo == True:
                if ref[i].isSmall == True:
                    ref_homo_small += 1
                    if found[i]:
                        homo_small += 1
                elif ref[i].isSmall == False:
                    ref_homo_large += 1
                    if found[i]:
                        homo_large += 1
            elif ref[i].isHomo == False:
                if ref[i].isSmall == True:
                    ref_hetero_small += 1
                    if found[i]:
                        hetero_small += 1
                elif ref[i].isSmall == False:
                    ref_hetero_large += 1
                    if found[i]:
                        hetero_large += 1
        print('Homozygous Small: %.3f Large: %.3f' %
              (homo_small / ref_homo_small * 100, homo_large / ref_homo_large * 100))
        print('Heterozygous Small: %.3f Large: %.3f' %
              (hetero_small/ref_hetero_small*100, hetero_large/ref_hetero_large*100))

    if PRINT_SPLIT and delsFlag:
        tool_large = tool_small = ref_large = ref_small = num_large = num_small = 0
        for i in range(len(tool)):
            len_cur = tool[i].en - tool[i].st + 1
            if len_cur >= SPLIT_LARGE:
                tool_large += 1
            else:
                tool_small += 1
        for i in range(len(ref)):
            len_cur = ref[i].en - ref[i].st + 1
            if len_cur >= SPLIT_LARGE:
                ref_large += 1
                if found[i]:
                    num_large += 1
            else:
                ref_small += 1
                if found[i]:
                    num_small += 1
        split_pr_large = float(100 * num_large / tool_large)
        split_re_large = float(100 * num_large / ref_large)
        split_f_large = float(
            2 * split_pr_large * split_re_large / (split_pr_large + split_re_large))
        split_pr_small = float(100 * num_small / tool_small)
        split_re_small = float(100 * num_small / ref_small)
        split_f_small = float(
            2 * split_pr_small * split_re_small / (split_pr_small + split_re_small))
        print('Small P: %.3f R: %.3f F: %.3f' %
              (split_pr_small, split_re_small, split_f_small))
        print('Large P: %.3f R: %.3f F: %.3f' %
              (split_pr_large, split_re_large, split_f_large))
    if GT_METRIC:
        print('GT: %.3f' % (gt_true/(gt_true + gt_false)))
    print('-' * 27)

    return dis, found, sup, le


def bpErrorPlot(bpe, toolLabels):
    fig, ax = plt.subplots()
    ax.set_ylabel('Breakpoint error (bp)')
    ax.boxplot(bpe, showfliers=False)
    ax.set_xticklabels(toolLabels)
    ax.set_ylim([-1, 6])


def supPlot(sup):
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    bins_sz = [2*x for x in range(25)]

    ax1.hist(sup[0], bins=bins_sz, ec='black')
    ax1.set_ylim([0, 200])
    ax1.set_ylabel('Number of deletions')
    ax1.set_xlabel('Breakpoint Support')
    ax1.set_title('HyINDEL')
    ax1.axvline(median(sup[0]), linestyle='dashed',
                color='k', linewidth=1)

    ax2.hist(sup[1], bins=bins_sz, ec='black')
    ax2.set_ylim([0, 200])
    ax2.set_ylabel('Number of deletions')
    ax2.set_xlabel('Breakpoint Support')
    ax2.set_title('Lumpy')
    ax2.axvline(median(sup[1]), linestyle='dashed',
                color='k', linewidth=1)

    ax3.hist(sup[2], bins=bins_sz, ec='black')
    ax3.set_ylim([0, 200])
    ax3.set_ylabel('Number of deletions')
    ax3.set_xlabel('Breakpoint Support')
    ax3.set_title('TIDDIT')
    ax3.axvline(median(sup[2]), linestyle='dashed',
                color='k', linewidth=1)

    ax4.hist(sup[3], bins=bins_sz, ec='black')
    ax4.set_ylim([0, 200])
    ax4.set_ylabel('Number of deletions')
    ax4.set_xlabel('Breakpoint Support')
    ax4.set_title('SoftSV')
    ax4.axvline(median(sup[3]), linestyle='dashed',
                color='k', linewidth=1)

    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    # fig.savefig('filename.eps', format='eps')


def lenErrorPlot(le):
    fig, ax = plt.subplots()
    ax.set_ylabel('Insertion Length error (bp)')
    ax.boxplot(le, showfliers=False)
    ax.set_xlabel('HyINDEL')


def cmp_common_ins(f1, f2):
    print(len(f1), len(f2))
    co = 0
    for i in range(len(f1)):
        ok = False
        for j in range(len(f2)):
            # if len(f2[j].seq) >= 50:
            if abs(f1[i].pos - f2[j].pos) <= 10:
                ok = True
        if ok:
            co += 1
    print(co)


def sim_stats(dels_ref, ins_ref):
    dels_homo_co = 0
    dels_hetero_co = 0
    for x in dels_ref:
        if x.isHomo:
            dels_homo_co += 1
        else:
            dels_hetero_co += 1
    print('Deletions Total:', len(dels_ref), '\tHomozygous:',
          dels_homo_co, '\tHeterozygous:', dels_hetero_co)
    ins_homo_co = 0
    ins_hetero_co = 0
    for x in ins_ref:
        if x.isHomo:
            ins_homo_co += 1
        else:
            ins_hetero_co += 1
    print('Insertions Total:', len(ins_ref), '\tHomozygous:',
          ins_homo_co, '\tHeterozygous:', ins_hetero_co)


def simulations():
    dels_ref = readSim('/Users/alok/IIIT/Simulations/Input/2.txt',
                       '/Users/alok/IIIT/Simulations/Input/1.txt', True)
    ins_ref = readSim('/Users/alok/IIIT/Simulations/Input/2.txt',
                      '/Users/alok/IIIT/Simulations/Input/1.txt', False)
    sim_stats(dels_ref, ins_ref)

    dels_hyindel_10x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/10x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_10x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/10x/hyindel/output.vcf', delsFlag=False)
    dels_hyindel_20x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/20x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_20x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/20x/hyindel/output.vcf', delsFlag=False)
    dels_hyindel_30x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/30x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_30x = readHyINDEL(
        '/Users/alok/IIIT/Simulations/Output/30x/hyindel/output.vcf', delsFlag=False)
    lumpy_10x = readLumpy(
        '/Users/alok/IIIT/Simulations/Output/10x/lumpy/lumpy_sim_10x.vcf')
    lumpy_20x = readLumpy(
        '/Users/alok/IIIT/Simulations/Output/20x/lumpy/lumpy_sim_20x.vcf')
    lumpy_30x = readLumpy(
        '/Users/alok/IIIT/Simulations/Output/30x/lumpy/lumpy_sim_30x.vcf')
    tiddit_10x = readTiddit(
        '/Users/alok/IIIT/Simulations/Output/10x/tiddit/tiddit_sim_10x.vcf')
    tiddit_20x = readTiddit(
        '/Users/alok/IIIT/Simulations/Output/20x/tiddit/tiddit_sim_20x.vcf')
    tiddit_30x = readTiddit(
        '/Users/alok/IIIT/Simulations/Output/30x/tiddit/tiddit_sim_30x.vcf')
    softsv_10x = readSoftsv('/Users/alok/IIIT/Simulations/Output/10x/softsv/deletions_small.txt',
                            '/Users/alok/IIIT/Simulations/Output/10x/softsv/deletions.txt')
    softsv_20x = readSoftsv('/Users/alok/IIIT/Simulations/Output/20x/softsv/deletions_small.txt',
                            '/Users/alok/IIIT/Simulations/Output/20x/softsv/deletions.txt')
    softsv_30x = readSoftsv('/Users/alok/IIIT/Simulations/Output/30x/softsv/deletions_small.txt',
                            '/Users/alok/IIIT/Simulations/Output/30x/softsv/deletions.txt')
    pamir_10x = readPamir(
        '/Users/alok/IIIT/Simulations/Output/10x/pamir_mrsfast_best/insertions_setcover.vcf')
    pamir_20x = readPamir(
        '/Users/alok/IIIT/Simulations/Output/20x/pamir_mrsfast_best/insertions_setcover.vcf')
    pamir_30x = readPamir(
        '/Users/alok/IIIT/Simulations/Output/30x/pamir_mrsfast_best/insertions_setcover.vcf')
    popins_10x = readPopins(
        '/Users/alok/IIIT/Simulations/Output/10x/popins/insertions.vcf')
    popins_20x = readPopins(
        '/Users/alok/IIIT/Simulations/Output/20x/popins/insertions.vcf')
    popins_30x = readPopins(
        '/Users/alok/IIIT/Simulations/Output/30x/popins/insertions.vcf')

    # Deletions
    compare(dels_hyindel_10x, dels_ref,
            'SIM-DELS-HYINDEL-10x', True, checkHomoAndSmall=True)
    compare(dels_hyindel_20x, dels_ref,
            'SIM-DELS-HYINDEL-20x', True, checkHomoAndSmall=True)
    b_hyindel_del, f, s_hyindel_del, l = compare(dels_hyindel_30x, dels_ref,
                                                 'SIM-DELS-HYINDEL-30x', True, checkHomoAndSmall=True)
    compare(lumpy_10x, dels_ref, 'SIM-DELS-LUMPY-10x', True)
    compare(lumpy_20x, dels_ref, 'SIM-DELS-LUMPY-20x', True)
    b_lumpy, f, s_lumpy, l = compare(
        lumpy_30x, dels_ref, 'SIM-DELS-LUMPY-30x', True)
    compare(tiddit_10x, dels_ref, 'SIM-DELS-TIDDIT-10x', True)
    compare(tiddit_20x, dels_ref, 'SIM-DELS-TIDDIT-20x', True)
    b_tiddit, f, s_tiddit, l = compare(
        tiddit_30x, dels_ref, 'SIM-DELS-TIDDIT-30x', True)
    compare(softsv_10x, dels_ref, 'SIM-DELS-SOFTSV-10x', True)
    compare(softsv_20x, dels_ref, 'SIM-DELS-SOFTSV-20x', True)
    b_softsv, f, s_softsv, l = compare(
        softsv_30x, dels_ref, 'SIM-DELS-SOFTSV-30x', True)

    # Insertions
    compare(ins_hyindel_10x, ins_ref,
            'SIM-INS-HYINDEL-10x', False, checkHomoAndSmall=True)
    compare(ins_hyindel_20x, ins_ref,
            'SIM-INS-HYINDEL-20x', False, checkHomoAndSmall=True)
    b_hyindel_ins, f, s, l_hyindel_ins = compare(ins_hyindel_30x, ins_ref,
                                                 'SIM-INS-HYINDEL-30x', False, checkHomoAndSmall=True)
    compare(pamir_10x, ins_ref, 'SIM-INS-PAMIR-10x', False)
    compare(pamir_20x, ins_ref, 'SIM-INS-PAMIR-20x', False)
    b_pamir, f, s, l = compare(pamir_30x, ins_ref, 'SIM-INS-PAMIR-30x', False)
    compare(popins_10x, ins_ref, 'SIM-INS-POPINS-10x', False)
    compare(popins_20x, ins_ref, 'SIM-INS-POPINS-20x', False)
    b_popins, f, s, l = compare(
        popins_30x, ins_ref, 'SIM-INS-POPINS-30x', False)

    # Plot breakpoint error
    if BREAKPOINT_PLOT:
        bpErrorPlot([b_hyindel_del, b_lumpy, b_tiddit, b_softsv],
                    ['HyINDEL', 'Lumpy', 'TIDDIT', 'SoftSV'])
        bpErrorPlot([b_hyindel_ins, b_pamir, b_popins],
                    ['HyINDEL', 'Pamir', 'Popins'])
        lenErrorPlot(l_hyindel_ins)
        plt.show()

    # Plot breakpoint support
    if SUPPORT_PLOT:
        print('Median breakpoint support (HyINDEL):', median(s_hyindel_del))
        print('Median breakpoint support (Lumpy):', median(s_lumpy))
        print('Median breakpoint support (Tiddit):', median(s_tiddit))
        print('Median breakpoint support (Softsv):', median(s_softsv))
        supPlot([s_hyindel_del, s_lumpy, s_tiddit, s_softsv])
        plt.show()


def platinum():
    dels_ref_svclassify = readSvclassify(
        '/Users/alok/IIIT/GS/Personalis_1000_Genomes_deduplicated_deletions.bed', True)
    ins_ref_svclassify = readSvclassify(
        '/Users/alok/IIIT/GS/Spiral_Genetics_insertions.bed', False)
    dels_ref_merged = readMergedBenchmark(
        '/Users/alok/IIIT/GS/NA12878_DGV-2016_LR-assembly.vcf', True)
    ins_ref_merged = readMergedBenchmark(
        '/Users/alok/IIIT/GS/NA12878_DGV-2016_LR-assembly.vcf', False)
    ins_ref_dgv_1kgp, ins_mob, ins = readDgvNovelIns(
        '/Users/alok/IIIT/GS/GRCh37_hg19_variants_2016-05-15.txt')

    dels_hyindel = readHyINDEL(
        '/Users/alok/IIIT/Platinum/Output/HyINDEL/output.vcf', delsFlag=True)
    ins_hyindel = readHyINDEL(
        '/Users/alok/IIIT/Platinum/Output/HyINDEL/output.vcf', delsFlag=False)
    lumpy = readLumpy(
        '/Users/alok/IIIT/Platinum/Output/lumpy/NA12878_lumpy.vcf')
    tiddit = readTiddit(
        '/Users/alok/IIIT/Platinum/Output/tiddit/NA12878_tiddit.vcf')
    softsv = readSoftsv('/Users/alok/IIIT/Platinum/Output/softsv/deletions_small.txt',
                        '/Users/alok/IIIT/Platinum/Output/softsv/deletions.txt')
    popins = readPopins(
        '/Users/alok/IIIT/Platinum/Output/popins/popins_platinum_insertions.vcf')
    pamir = readPamir(
        '/Users/alok/IIIT/Platinum/Output/pamir/platinum_pamir_insertions_setcover.vcf')

    # Deletions svclassify
    compare(dels_hyindel, dels_ref_svclassify,
            'PLATINUM-DELS-HYINDEL-SVCLASSIFY', True, checkSmall=True)
    compare(lumpy, dels_ref_svclassify,
            'PLATINUM-DELS-LUMPY-SVCLASSIFY', True, checkSmall=True)
    compare(tiddit, dels_ref_svclassify,
            'PLATINUM-DELS-TIDDIT-SVCLASSIFY', True, checkSmall=True)
    compare(softsv, dels_ref_svclassify,
            'PLATINUM-DELS-SOFTSV-SVCLASSIFY', True, checkSmall=True)

    # Deletions merged
    compare(dels_hyindel, dels_ref_merged,
            'PLATINUM-DELS-HYINDEL-MERGED', True, checkSmall=True)
    compare(lumpy, dels_ref_merged,
            'PLATINUM-DELS-LUMPY-MERGED', True, checkSmall=True)
    compare(tiddit, dels_ref_merged,
            'PLATINUM-DELS-TIDDIT-MERGED', True, checkSmall=True)
    compare(softsv, dels_ref_merged,
            'PLATINUM-DELS-SOFTSV-MERGED', True, checkSmall=True)

    # Insertions dgv novel sequence
    compare(ins_hyindel, ins_ref_dgv_1kgp,
            'PLATINUM-INS-HYINDEL-DGVNOVELINS', False)
    compare(pamir, ins_ref_dgv_1kgp, 'PLATINUM-INS-PAMIR-DGVNOVELINS', False)
    compare(popins, ins_ref_dgv_1kgp, 'PLATINUM-INS-POPINS-DGVNOVELINS', False)

    # Insertions svclassify
    compare(ins_hyindel, ins_ref_svclassify,
            'PLATINUM-INS-HYINDEL-SVCLASSIFY', False)
    compare(pamir, ins_ref_svclassify, 'PLATINUM-INS-PAMIR-SVCLASSIFY', False)
    compare(popins, ins_ref_svclassify,
            'PLATINUM-INS-POPINS-SVCLASSIFY', False)

    # Insertions merged
    compare(ins_hyindel, ins_ref_merged, 'PLATINUM-INS-HYINDEL-MERGED', False)
    compare(pamir, ins_ref_merged, 'PLATINUM-INS-PAMIR-MERGED', False)
    compare(popins, ins_ref_merged, 'PLATINUM-INS-POPINS-MERGED', False)


def main():
    simulations()
    # platinum()
    return


if __name__ == '__main__':
    main()
