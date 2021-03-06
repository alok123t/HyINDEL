# pylint: disable=unused-variable

from statistics import mean, median
import matplotlib.pyplot as plt

SPLIT_LARGE = 500
# Remove events < 50 bp in tools (True)
FLAG_50 = True
# Plot fscore (True)
FSCORE_PLOT = False
# Plot breakpoint error plot (True)
BREAKPOINT_PLOT = True
# Plot support (True)
SUPPORT_PLOT = True
# Plot size distributions (True)
SIZE_DISTR_PLOT = False
# Maximum error in breakpoint for insertions
DIST_INS = 10
# Reciprocal overlap
RO = 0.5


class Deletion:
    def __init__(self, t_chr, t_st, t_en, t_iD=-1, t_pe=-1, t_sc=-1, t_sr=-1, t_isHomo="N/A", t_isSmall="N/A", t_isImprecise="N/A"):
        self.chr = modChr(t_chr)
        self.st = t_st
        self.en = t_en
        self.iD = t_iD
        self.pe = t_pe
        self.sc = t_sc
        self.sr = t_sr
        self.isHomo = t_isHomo
        self.isSmall = t_isSmall
        self.isImprecise = t_isImprecise

    def __str__(self):
        return self.chr+'\t' + str(self.st)+'\t' + str(self.en)+'\t' + str(self.isImprecise)


class Insertion:
    def __init__(self, t_chr, t_pos, t_seq='', t_seqLen=-1, t_refPos=-1, t_iD=-1, t_sc=-1, t_sup=-1, t_isHomo="N/A", t_isSmall="N/A", t_isImprecise="N/A"):
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
        self.isImprecise = t_isImprecise

    def __str__(self):
        if len(self.seq) == 0:
            return self.chr + '\t' + str(self.pos) + str(self.isImprecise)
        else:
            return self.chr + '\t' + str(self.pos) + '\t' + str(len(self.seq)) + str(self.isImprecise)


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
            isImprecise = False
            if 'IMPRECISE' in line:
                isImprecise = True
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

                st = int(l[1])
                en = int(info_col[2].split('=')[1])
                sz = en - st + 1
                isSmall = False
                if sz <= 500:
                    isSmall = True
                ret.append(Deletion(t_chr=l[0], t_st=st, t_en=en, t_pe=sup_pe, t_sc=sup_sc,
                                    t_sr=sup_sr, t_isHomo=isHomo, t_isSmall=isSmall, t_isImprecise=isImprecise))
            else:
                if l[4] != '<INS>':
                    continue
                if 'IMPRECISE' in l[7] or 'SEQUP' in l[7] or 'SEQDOWN' in l[7]:
                    ret.append(Insertion(t_chr=l[0], t_pos=int(
                        l[1]), t_isImprecise=isImprecise))
                else:
                    info_col = l[7].split(';')
                    seq = info_col[2].split('=')[1]
                    seqLen = len(seq)
                    isSmall = False
                    if seqLen <= 500:
                        isSmall = True
                    ret.append(Insertion(t_chr=l[0], t_pos=int(
                        l[1]), t_seq=seq, t_seqLen=seqLen, t_isSmall=isSmall, t_isImprecise=isImprecise))

    return ret


def readLumpy(fName):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split()
            isImprecise = False
            if 'IMPRECISE' in line:
                isImprecise = True
            if l[4] != '<DEL>':
                continue
            info_col = l[7].split(';')
            format_col = l[9].split(':')
            sup_pe = int(format_col[2])
            sup_sr = int(format_col[3])
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(info_col[3].split(
                '=')[1]), t_pe=sup_pe, t_sr=sup_sr, t_isImprecise=isImprecise))

    return ret


def readTiddit(fName):
    ret = []
    with open(fName, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            l = line.split()
            isImprecise = False
            if 'ORSR=0,0' in line:
                isImprecise = True
            if l[4] != '<DEL>':
                continue
            # if l[6] != 'PASS':
            #     continue
            info_col = l[7].split(';')
            format_col = l[9].split(':')
            sup_pe = int(format_col[2])
            sup_sr = int(format_col[3])
            quala = int(info_col[-2].split('=')[1])
            qualb = int(info_col[-1].split('=')[1])
            if quala <= 55 or qualb <= 55:
                isImprecise = True
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                info_col[3].split('=')[1]), t_pe=sup_pe, t_sr=sup_sr, t_isImprecise=isImprecise))

    return ret


def readSoftsv(fNameSmall, fNameLarge):
    ret = []
    with open(fNameLarge, 'r') as fLarge, open(fNameSmall, 'r') as fSmall:
        for line in fSmall:
            l = line.split()
            if l[0] == 'Chromosome':
                continue
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                l[2]), t_pe=int(l[4]), t_sc=int(l[5]), t_isImprecise=False))

        for line in fLarge:
            l = line.split()
            if l[0] == 'Chromosome':
                continue
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
                l[2]), t_pe=int(l[4]), t_sc=int(l[5]), t_isImprecise=False))
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
            l[1]), t_seq=seq, t_sup=sup, t_isImprecise=False))

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
            ret.append(Insertion(t_chr=l[0], t_pos=int(
                l[1]), t_isImprecise=False))

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


def readDgv(fName, delsFlag):
    ret = []
    ret_mob = []
    ret_ins = []
    ret_loss = []
    ret_dels = []
    with open(fName, 'r') as f:
        for line in f:
            l = line.split('\t')
            if l[0] == 'variantaccession':
                continue
            if 'NA12878' not in line:
                continue
            if delsFlag:
                st = int(l[2])
                en = int(l[3])
                sz = en - st + 1
                isSmall = False
                if sz <= 500:
                    isSmall = True
                if l[5] == 'deletion':
                    if l[6] == '1000_Genomes_Consortium_Phase_1':
                        ret_dels.append(
                            Deletion(t_chr=l[1], t_st=st, t_en=en, t_isSmall=isSmall))
                if l[5] == 'loss':
                    if l[6] == '1000_Genomes_Consortium_Phase_3':
                        ret_loss.append(
                            Deletion(t_chr=l[1], t_st=st, t_en=en, t_isSmall=isSmall))
            else:
                if l[6] != '1000_Genomes_Consortium_Pilot_Project':
                    continue
                if l[5] == 'novel sequence insertion':
                    ret.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
                if l[5] == 'mobile element insertion':
                    ret_mob.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
                if l[5] == 'insertion':
                    ret_ins.append(Insertion(t_chr=l[1], t_pos=int(l[2])))
    if delsFlag:
        return ret_dels, ret_loss
    else:
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
                isSmall = False
                if le <= 500:
                    isSmall = True
                ret.append(Insertion(t_chr=l[0], t_pos=int(
                    l[1]), t_seqLen=le, t_isSmall=isSmall))

    f.close()
    return ret


"""
returns [[breakpoint error], [found], [support], [length error]]
[breakpoint error]: only for detected variants
[found]: True/False if present in benchmark
[support]: sc/sr support only for detected variants (if present)
[length error]: only for precise insertions detected variants
[fScore, homo_small_recall, homo_large_recall, hetero_small_recall, hetero_large_recall]
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

    found = [-1] * len(tool)
    print(len(tool))

    dis = []
    sup = []
    le = []
    gt_true = 0
    gt_false = 0
    for i in range(len(tool)):
        for j in range(len(ref)):
            if delsFlag:
                if checkDel(tool[i], ref[j]):
                    found[i] = j
                    if tool[i].isImprecise == False:
                        dis.append(abs(tool[i].en - ref[j].en) +
                                   abs(tool[i].st - ref[j].st))
                        if dis[-1] >= 20:
                            print(tool[i], ref[j])
                    exact_sup = tool[i].sc + tool[i].sr
                    if exact_sup > 0:
                        sup.append(exact_sup)
            else:
                if checkIns(tool[i], ref[j]):
                    found[i] = j
                    if tool[i].isImprecise == False:
                        dis.append(abs(tool[i].pos - ref[j].pos))
                    if tool[i].seqLen != -1:
                        le.append(ref[j].seqLen - tool[i].seqLen)

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = len(set(found))

    precision = float(100 * num_true / num_tool)
    recall = float(100 * num_true / num_ref)
    fScore = float(2 * precision*recall/(precision+recall))
    print('P: %.3f R: %.3f F: %.3f' % (precision, recall, fScore))
    retVals = []
    retVals.append(fScore)

    if not delsFlag:
        print('Found:', num_true)

    # For real data
    if checkSmall:
        ref_small = 0
        ref_large = 0
        small = []
        large = []
        for i in range(len(ref)):
            if ref[i].isSmall == True:
                ref_small += 1
            elif ref[i].isSmall == False:
                ref_large += 1
        for i in range(len(found)):
            jIdx = found[i]
            if jIdx != -1:
                if ref[jIdx].isSmall == True:
                    small.append(jIdx)
                elif ref[jIdx].isSmall == False:
                    large.append(jIdx)
        small_co = len(set(small))
        large_co = len(set(large))

        print('Small: %.3f Large: %.3f' %
              (small_co / ref_small * 100, large_co / ref_large * 100))

    # For simulated data
    if checkHomoAndSmall:
        ref_homo_small = 0
        ref_homo_large = 0
        ref_hetero_small = 0
        ref_hetero_large = 0
        for i in range(len(ref)):
            if ref[i].isHomo == True:
                if ref[i].isSmall == True:
                    ref_homo_small += 1
                elif ref[i].isSmall == False:
                    ref_homo_large += 1
            elif ref[i].isHomo == False:
                if ref[i].isSmall == True:
                    ref_hetero_small += 1
                elif ref[i].isSmall == False:
                    ref_hetero_large += 1

        homo_small = []
        homo_large = []
        hetero_small = []
        hetero_large = []
        for i in range(len(found)):
            jIdx = found[i]
            if jIdx != -1:
                if ref[jIdx].isHomo == True:
                    if ref[jIdx].isSmall == True:
                        homo_small.append(jIdx)
                    elif ref[jIdx].isSmall == False:
                        homo_large.append(jIdx)
                elif ref[jIdx].isHomo == False:
                    if ref[jIdx].isSmall == True:
                        hetero_small.append(jIdx)
                    elif ref[jIdx].isSmall == False:
                        hetero_large.append(jIdx)

        homo_small_co = len(set(homo_small))
        homo_large_co = len(set(homo_large))
        hetero_small_co = len(set(hetero_small))
        hetero_large_co = len(set(hetero_large))
        print('Homozygous Small: %.3f Large: %.3f' %
              (homo_small_co / ref_homo_small * 100, homo_large_co / ref_homo_large * 100))
        print('Heterozygous Small: %.3f Large: %.3f' %
              (hetero_small_co / ref_hetero_small*100, hetero_large_co / ref_hetero_large*100))
        retVals.append(homo_small_co / ref_homo_small * 100)
        retVals.append(homo_large_co / ref_homo_large * 100)
        retVals.append(hetero_small_co / ref_hetero_small*100)
        retVals.append(hetero_large_co / ref_hetero_large*100)

    print('-' * 27)

    return dis, found, sup, le, retVals


def recallPlot(f):
    fig, ax = plt.subplots()
    x = [i for i in range(3)]
    ax.set_xticks(x)
    ax.set_xticklabels(['10x', '20x', '30x'])
    ax.set_xlabel('Coverage')
    ax.set_ylim([0, 100])
    ax.set_ylabel('Recall')

    homo_small = [f[0][1], f[1][1], f[2][1]]
    homo_large = [f[0][2], f[1][2], f[2][2]]
    hetero_small = [f[0][3], f[1][3], f[2][3]]
    hetero_large = [f[0][4], f[1][4], f[2][4]]

    ax.plot(x, homo_small, label='Homozygous small', color='c', marker='o')
    ax.plot(x, homo_large, label='Homozygous large', color='r', marker='v')
    ax.plot(x, hetero_small, label='Heterozygous small', color='b', marker='s')
    ax.plot(x, hetero_large, label='Heterozygous large', color='g', marker='D')

    ax.legend(loc='lower left')


def fScorePlot(f, delsFlag):
    fig, ax = plt.subplots()
    x = [i for i in range(3)]
    ax.set_xticks(x)
    ax.set_xticklabels(['10x', '20x', '30x'])
    ax.set_xlabel('Coverage')
    ax.set_ylim([0, 100])
    ax.set_ylabel('F-score')

    if delsFlag:
        ax.plot(x, f[0], label='HyINDEL', color='c', marker='o')
        ax.plot(x, f[1], label='Lumpy', color='r', marker='v')
        ax.plot(x, f[2], label='TIDDIT', color='b', marker='s')
        ax.plot(x, f[3], label='SoftSV', color='g', marker='D')
    else:
        ax.plot(x, f[0], label='HyINDEL', color='c', marker='o')
        ax.plot(x, f[1], label='Pamir', color='r', marker='v')
        ax.plot(x, f[2], label='Popins', color='b', marker='s')

    ax.legend(loc='lower left')


def bpErrorPlot(bpe, toolLabels):
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['font.size'] = 14
    fig, ax = plt.subplots()
    # Deletions
    if len(bpe) == 4:
        ax.set_title('(a)', fontweight='bold', position=(0.9, 0.9))
        ax.set_ylim(ymin=-2, ymax=12)
    # Insertions
    elif len(bpe) == 3:
        ax.set_title('(b)', fontweight='bold', position=(0.9, 0.9))
        ax.set_ylim(ymin=-2, ymax=10)
    ax.boxplot(bpe)
    ax.set_xticklabels(toolLabels, fontweight='bold')
    ax.set_ylabel('Breakpoint error (bp)', fontweight='bold')


def supPlot(sup, toolLabels):
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['font.size'] = 12
    fig, ax = plt.subplots(
        1, 4, sharey=True, figsize=(8, 4))
    fig.text(0.5, 0.04, 'Breakpoint support', ha='center')

    bins_sz = [2 * x for x in range(25)]

    for i in range(4):
        if i == 0:
            ax[i].set_ylabel('Number of deletions')
        ax[i].hist(sup[i], bins=bins_sz, ec='black')
        ax[i].set_ylim([0, 200])
        ax[i].set_title(toolLabels[i], fontweight='bold', position=(0.8, 0.85))
        ax[i].axvline(median(sup[i]), linestyle='dashed',
                      color='k', linewidth=1)

    plt.tight_layout()
    plt.subplots_adjust(left=0.10, bottom=0.14, right=0.98)


def lenErrorPlot(le):
    fig, ax = plt.subplots()
    ax.boxplot(le)
    ax.set_xlabel('HyINDEL')
    ax.set_ylabel('Insertion Length error (bp)')
    # ax.set_ylim([-10, 100])


def sizeDistrPlot(f, delsFlag):
    if delsFlag:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ar = []
        for x in f:
            ar.append(x.en - x.st + 1)

        bins_sz1 = [20 * x for x in range(50)]
        ax1.set_xlabel('Small deletions size (bp)')
        ax1.set_ylabel('Frequency')
        ax1.hist(ar, bins=bins_sz1, ec='black')

        bins_sz2 = [1000 + (200 * x) for x in range(50)]
        ax2.set_xlabel('Large deletions size (bp)')
        ax2.set_ylabel('Frequency')
        ax2.hist(ar, bins=bins_sz2, ec='black')
    else:
        fig, ax = plt.subplots()
        ar = []
        for x in f:
            if x.seqLen != -1:
                ar.append(x.seqLen)

        bins_sz = [20 * x for x in range(100)]
        ax.set_xlabel('Insertions size (bp)')
        ax.set_ylabel('Frequency')
        ax.hist(ar, bins=bins_sz, ec='black')


def cmp_common_ins2(f1, f2, ex):
    print(len(f1), len(f2))
    co = 0
    for i in range(len(f1)):
        ok = False
        rem = False
        for j in range(len(f2)):
            if f1[i].chr != f2[j].chr:
                continue
            if abs(f1[i].pos - f2[j].pos) <= 10:
                for k in range(len(ex)):
                    if abs(f1[i].pos - ex[k].pos) <= 10:
                        rem = True
                        break
                ok = True
                break
        if ok and not rem:
            co += 1
            # print(f1[i])
    print('Common:', co)


def cmp_common_ins3(f1, f2, f3):
    print('HyINDEL:', len(f1), 'Pamir:', len(f2), 'Popins:', len(f3))
    l = []
    for i in range(len(f1)):
        ok = False
        for j in range(len(f2)):
            if abs(f1[i].pos - f2[j].pos) <= 10:
                ok = True
                break
        if ok:
            l.append(f1[i])
    co = 0
    ret = []
    for i in range(len(l)):
        ok = False
        for j in range(len(f3)):
            if abs(l[i].pos - f3[j].pos) <= 10:
                ret.append(l[i])
                ok = True
                break
        if ok:
            co += 1
    print('Common3:', co)
    return ret


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


def real_stats(f, extra=False, delsFlag=False):
    small_co = 0
    large_co = 0
    imprecise_co = 0
    ar = []
    for x in f:
        if extra:
            if delsFlag:
                ar.append(x.en-x.st+1)
            else:
                if x.seqLen != -1:
                    ar.append(x.seqLen)

        if x.isSmall == True:
            small_co += 1
        elif x.isSmall == False:
            large_co += 1
        else:
            imprecise_co += 1
    print('Small:', small_co, 'Large:', large_co,
          'Imprecise:', imprecise_co)
    if extra:
        ar.sort()
        print('Min:', ar[0], 'Max:', ar[-1], 'Mean:',
              mean(ar), 'Median:', median(ar))


def ref_stats(dels_ref_svclassify, ins_ref_svclassify, dels_ref_merged, ins_ref_merged, ins_ref_dgv_1kgp, dels_ref_dgv_1kgp_loss):
    print('svclassify-dels')
    real_stats(dels_ref_svclassify)
    print('Total:', len(dels_ref_svclassify))
    print('-' * 27)
    print('svclassify-ins')
    print('Total:', len(ins_ref_svclassify))
    print('-' * 27)
    print('merged-dels')
    real_stats(dels_ref_merged)
    print('Total:', len(dels_ref_merged))
    print('-' * 27)
    print('merged-ins')
    real_stats(ins_ref_merged)
    print('Total:', len(ins_ref_merged))
    print('-' * 27)
    print('dgv-novelins')
    print('Total:', len(ins_ref_dgv_1kgp))
    print('-' * 27)
    print('dgv-dels')
    real_stats(dels_ref_dgv_1kgp_loss)
    print('Total:', len(dels_ref_dgv_1kgp_loss))
    print('-' * 27)


def tool_stats(dels_hyindel, ins_hyindel, pamir, popins):
    print('hyindel-dels:', len(dels_hyindel))
    real_stats(dels_hyindel, extra=True, delsFlag=True)
    if SIZE_DISTR_PLOT:
        sizeDistrPlot(dels_hyindel, delsFlag=True)
        plt.subplots_adjust(wspace=0.3)
        sizeDistrPlot(ins_hyindel, delsFlag=False)
        plt.show()
    print('hyindel-ins:', len(ins_hyindel))
    real_stats(ins_hyindel, extra=True, delsFlag=False)
    print('pamir:', len(pamir))
    print('popins:', len(popins))
    print('-' * 27)

    ex = cmp_common_ins3(ins_hyindel, pamir, popins)
    print('-' * 27)
    ex = []  # do not exclude any
    cmp_common_ins2(ins_hyindel, pamir, ex)
    print('-' * 27)
    cmp_common_ins2(ins_hyindel, popins, ex)
    print('-' * 27)
    cmp_common_ins2(pamir, popins, ex)
    print('-' * 27)


def simulations():
    dels_ref = readSim('/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Input/2.txt',
                       '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Input/1.txt', True)
    ins_ref = readSim('/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Input/2.txt',
                      '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Input/1.txt', False)
    # sim_stats(dels_ref, ins_ref)

    dels_hyindel_10x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_10x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/hyindel/output.vcf', delsFlag=False)
    dels_hyindel_20x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_20x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/hyindel/output.vcf', delsFlag=False)
    dels_hyindel_30x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/hyindel/output.vcf', delsFlag=True)
    ins_hyindel_30x = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/hyindel/output.vcf', delsFlag=False)
    lumpy_10x = readLumpy(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/lumpy/lumpy_sim_10x.vcf')
    lumpy_20x = readLumpy(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/lumpy/lumpy_sim_20x.vcf')
    lumpy_30x = readLumpy(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/lumpy/lumpy_sim_30x.vcf')
    tiddit_10x = readTiddit(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/tiddit/tiddit_sim_10x.vcf')
    tiddit_20x = readTiddit(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/tiddit/tiddit_sim_20x.vcf')
    tiddit_30x = readTiddit(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/tiddit/tiddit_sim_30x.vcf')
    softsv_10x = readSoftsv('/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/softsv/deletions_small.txt',
                            '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/softsv/deletions.txt')
    softsv_20x = readSoftsv('/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/softsv/deletions_small.txt',
                            '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/softsv/deletions.txt')
    softsv_30x = readSoftsv('/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/softsv/deletions_small.txt',
                            '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/softsv/deletions.txt')
    pamir_10x = readPamir(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/pamir_mrsfast_best/insertions_setcover.vcf')
    pamir_20x = readPamir(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/pamir_mrsfast_best/insertions_setcover.vcf')
    pamir_30x = readPamir(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/pamir_mrsfast_best/insertions_setcover.vcf')
    popins_10x = readPopins(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/10x/popins/insertions.vcf')
    popins_20x = readPopins(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/20x/popins/insertions.vcf')
    popins_30x = readPopins(
        '/Volumes/GoogleDrive/My Drive/IIIT/Simulations/Output/30x/popins/insertions.vcf')

    # Deletions
    out_dels_hyindel_10x = compare(dels_hyindel_10x, dels_ref,
                                   'SIM-DELS-HYINDEL-10x', True, checkHomoAndSmall=True)
    out_dels_hyindel_20x = compare(dels_hyindel_20x, dels_ref,
                                   'SIM-DELS-HYINDEL-20x', True, checkHomoAndSmall=True)
    out_dels_hyindel_30x = compare(dels_hyindel_30x, dels_ref,
                                   'SIM-DELS-HYINDEL-30x', True, checkHomoAndSmall=True)
    out_dels_lumpy_10x = compare(
        lumpy_10x, dels_ref, 'SIM-DELS-LUMPY-10x', True)
    out_dels_lumpy_20x = compare(
        lumpy_20x, dels_ref, 'SIM-DELS-LUMPY-20x', True)
    out_dels_lumpy_30x = compare(
        lumpy_30x, dels_ref, 'SIM-DELS-LUMPY-30x', True)
    out_dels_tiddit_10x = compare(
        tiddit_10x, dels_ref, 'SIM-DELS-TIDDIT-10x', True)
    out_dels_tiddit_20x = compare(
        tiddit_20x, dels_ref, 'SIM-DELS-TIDDIT-20x', True)
    out_dels_tiddit_30x = compare(
        tiddit_30x, dels_ref, 'SIM-DELS-TIDDIT-30x', True)
    out_dels_softsv_10x = compare(
        softsv_10x, dels_ref, 'SIM-DELS-SOFTSV-10x', True)
    out_dels_softsv_20x = compare(
        softsv_20x, dels_ref, 'SIM-DELS-SOFTSV-20x', True)
    out_dels_softsv_30x = compare(
        softsv_30x, dels_ref, 'SIM-DELS-SOFTSV-30x', True)

    # Insertions
    out_ins_hyindel_10x = compare(ins_hyindel_10x, ins_ref,
                                  'SIM-INS-HYINDEL-10x', False, checkHomoAndSmall=True)
    out_ins_hyindel_20x = compare(ins_hyindel_20x, ins_ref,
                                  'SIM-INS-HYINDEL-20x', False, checkHomoAndSmall=True)
    out_ins_hyindel_30x = compare(ins_hyindel_30x, ins_ref,
                                  'SIM-INS-HYINDEL-30x', False, checkHomoAndSmall=True)
    out_ins_pamir_10x = compare(
        pamir_10x, ins_ref, 'SIM-INS-PAMIR-10x', False)
    out_ins_pamir_20x = compare(
        pamir_20x, ins_ref, 'SIM-INS-PAMIR-20x', False)
    out_ins_pamir_30x = compare(
        pamir_30x, ins_ref, 'SIM-INS-PAMIR-30x', False)
    out_ins_popins_10x = compare(
        popins_10x, ins_ref, 'SIM-INS-POPINS-10x', False)
    out_ins_popins_20x = compare(
        popins_20x, ins_ref, 'SIM-INS-POPINS-20x', False)
    out_ins_popins_30x = compare(
        popins_30x, ins_ref, 'SIM-INS-POPINS-30x', False)

    # Plot fscore
    if FSCORE_PLOT:
        recallPlot([out_dels_hyindel_10x[4],
                    out_dels_hyindel_20x[4], out_dels_hyindel_30x[4]])
        recallPlot([out_ins_hyindel_10x[4],
                    out_ins_hyindel_20x[4], out_ins_hyindel_30x[4]])
        # fScorePlot([[f_hyindel_del_10x[0], f_hyindel_del_20x[0], f_hyindel_del_30x[0]], [f_lumpy_10x[0], f_lumpy_20x[0], f_lumpy_30x[0]], [
        #     f_tiddit_10x[0], f_tiddit_20x[0], f_tiddit_30x[0]], [f_softsv_10x[0], f_softsv_20x[0], f_softsv_30x[0]]], delsFlag=True)
        # fScorePlot([[f_hyindel_ins_10x[0], f_hyindel_ins_20x[0], f_hyindel_ins_30x[0]], [
        #            f_pamir_10x[0], f_pamir_20x[0], f_pamir_30x[0]], [f_popins_10x[0], f_popins_20x[0], f_popins_30x[0]]], delsFlag=False)

    # Plot breakpoint error
    if BREAKPOINT_PLOT:
        bpErrorPlot([out_dels_hyindel_30x[0], out_dels_lumpy_30x[0], out_dels_tiddit_30x[0], out_dels_softsv_30x[0]],
                    ['(i)', '(ii)', '(iii)', '(iv)'])
        bpErrorPlot([out_ins_hyindel_30x[0], out_ins_pamir_30x[0], out_ins_popins_30x[0]],
                    ['(i)', '(ii)', '(iii)'])
        # lenErrorPlot(l_hyindel_ins)

    # Plot breakpoint support
    if SUPPORT_PLOT:
        print('Median breakpoint support (HyINDEL):',
              median(out_dels_hyindel_30x[2]))
        print('Median breakpoint support (Lumpy):',
              median(out_dels_lumpy_30x[2]))
        print('Median breakpoint support (Tiddit):',
              median(out_dels_tiddit_30x[2]))
        print('Median breakpoint support (Softsv):',
              median(out_dels_softsv_30x[2]))
        supPlot([out_dels_hyindel_30x[2], out_dels_lumpy_30x[2], out_dels_tiddit_30x[2], out_dels_softsv_30x[2]],
                ['(i)', '(ii)', '(iii)', '(iv)'])

    if FSCORE_PLOT or BREAKPOINT_PLOT or SUPPORT_PLOT:
        plt.show()


def platinum():
    dels_ref_svclassify = readSvclassify(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/Personalis_1000_Genomes_deduplicated_deletions.bed', True)
    ins_ref_svclassify = readSvclassify(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/Spiral_Genetics_insertions.bed', False)
    dels_ref_merged = readMergedBenchmark(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/NA12878_DGV-2016_LR-assembly.vcf', True)
    ins_ref_merged = readMergedBenchmark(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/NA12878_DGV-2016_LR-assembly.vcf', False)
    dels_ref_dgv_1kgp, dels_ref_dgv_1kgp_loss = readDgv(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/GRCh37_hg19_variants_2016-05-15.txt', True)
    ins_ref_dgv_1kgp, ins_mob, ins = readDgv(
        '/Volumes/GoogleDrive/My Drive/IIIT/Benchmarks/GRCh37_hg19_variants_2016-05-15.txt', False)

    dels_hyindel = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/hyindel/output.vcf', delsFlag=True)
    ins_hyindel = readHyINDEL(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/hyindel/output.vcf', delsFlag=False)
    lumpy = readLumpy(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/lumpy/NA12878_lumpy.vcf')
    tiddit = readTiddit(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/tiddit/NA12878_tiddit.vcf')
    softsv = readSoftsv('/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/softsv/deletions_small.txt',
                        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/softsv/deletions.txt')
    popins = readPopins(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/popins/popins_platinum_insertions.vcf')
    pamir = readPamir(
        '/Volumes/GoogleDrive/My Drive/IIIT/Platinum/Output/pamir/platinum_pamir_insertions_setcover.vcf')

    ref_stats(dels_ref_svclassify, ins_ref_svclassify, dels_ref_merged,
              ins_ref_merged, ins_ref_dgv_1kgp, dels_ref_dgv_1kgp_loss)
    tool_stats(dels_hyindel, ins_hyindel, pamir, popins)

    # Deletions dgv
    compare(dels_hyindel, dels_ref_dgv_1kgp,
            'PLATINUM-DELS-HYINDEL-DGV', True)
    compare(lumpy, dels_ref_dgv_1kgp,
            'PLATINUM-DELS-LUMPY-DGV', True)
    compare(tiddit, dels_ref_dgv_1kgp,
            'PLATINUM-DELS-TIDDIT-DGV', True)
    compare(softsv, dels_ref_dgv_1kgp,
            'PLATINUM-DELS-SOFTSV-DGV', True)

    # Deletions dgv loss
    compare(dels_hyindel, dels_ref_dgv_1kgp_loss,
            'PLATINUM-DELS-HYINDEL-DGV-LOSS', True, checkSmall=True)
    compare(lumpy, dels_ref_dgv_1kgp_loss,
            'PLATINUM-DELS-LUMPY-DGV-LOSS', True, checkSmall=True)
    compare(tiddit, dels_ref_dgv_1kgp_loss,
            'PLATINUM-DELS-TIDDIT-DGV-LOSS', True, checkSmall=True)
    compare(softsv, dels_ref_dgv_1kgp_loss,
            'PLATINUM-DELS-SOFTSV-DGV-LOSS', True, checkSmall=True)

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
    plt.rcParams['axes.labelweight'] = 'bold'

    simulations()
    platinum()


if __name__ == '__main__':
    main()
