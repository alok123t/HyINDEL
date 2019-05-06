# pylint: disable=unused-variable

from statistics import median
import numpy as np
import matplotlib.pyplot as plt

# Remove events < 50 bp in tools
FLAG_50 = True
# Plot breakpoint error plot
BREAKPOINT_PLOT = False
# Plot support
SUPPORT_PLOT = False
# Print metrics seperately for large and small variants
PRINT_SPLIT = True
SPLIT_LARGE = 500
# Maximum error in breakpoint for insertions
DIST_INS = 10
# Reciprocal overlap
SIM_RO = 0.5
REAL_RO = 0.5
GT_METRIC = False

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


def checkDel(a, b, RO):
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


class Deletion:
    def __init__(self, t_chr, t_st, t_en, t_iD=-1, t_pe=-1, t_sc=-1, t_sr=-1, t_isHomo="N/A"):
        self.chr = modChr(t_chr)
        self.st = t_st
        self.en = t_en
        self.iD = t_iD
        self.pe = t_pe
        self.sc = t_sc
        self.sr = t_sr
        self.isHomo = t_isHomo

    def __str__(self):
        return self.chr+'\t' + str(self.st)+'\t' + str(self.en)


class Insertion:
    def __init__(self, t_chr, t_pos, t_seq='', t_seqLen=-1, t_refPos=-1, t_iD=-1, t_sc=-1, t_isHomo=True):
        self.chr = modChr(t_chr)
        self.pos = t_pos  # Position on reference
        self.seq = t_seq
        self.seqLen = t_seqLen
        self.refPos = t_refPos  # Position on sample
        self.iD = t_iD
        self.sc = t_sc
        self.isHomo = t_isHomo

    def __str__(self):
        return self.chr+'\t' + str(self.pos)


def readMinus(fName):
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
        info_col = l[7].split(';')
        format_col = l[9].split(':')
        sup_pe = int(format_col[2])
        sup_sr = int(format_col[3])
        ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
            info_col[3].split('=')[1]), t_pe=sup_pe, t_sr=sup_sr))

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
        info_col = l[7].split(';')
        format_col = l[9].split(':')
        sup_pe = int(format_col[2])
        sup_sr = int(format_col[3])
        ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
            info_col[3].split('=')[1]), t_pe=sup_pe, t_sr=sup_sr))

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
        ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
            l[2]), t_pe=int(l[4]), t_sc=int(l[5])))

    inp = open(fNameLarge)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'Chromosome':
            continue
        ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(
            l[2]), t_pe=int(l[4]), t_sc=int(l[5])))

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
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(l[2])))
        else:
            ret.append(Insertion(t_chr=l[0], t_pos=int(l[1])))

    return ret


def readSim(fName, delsFlag):
    inp = open(fName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for i in range(len(lines)):
        l = lines[i].split()

        isHomo = True
        if i >= 1000:
            isHomo = False
        if l[0] == 'DELL':
            if delsFlag:
                ret.append(Deletion(t_chr=l[1], t_st=int(
                    l[2]), t_en=int(l[3]), t_isHomo=isHomo))
        else:
            if not delsFlag:
                ret.append(Insertion(t_chr=l[1], t_pos=int(
                    l[2]), t_seq=l[4]), t_isHomo=isHomo)

    return ret


"""
returns [[breakpoint error], [found], [support]]
[breakpoint error]: only for detected variants
[found]: True/False if present in benchmark
[support]: sc/sr support only for detected variants (if present)
"""


def compare(toolUnfiltered, ref, toolName, delsFlag, RO=0.5):
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
    gt_true = 0
    gt_false = 0
    for i in range(len(tool)):
        for j in range(len(ref)):
            if delsFlag:
                if checkDel(tool[i], ref[j], RO):
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
                    print(dis[-1])

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = sum(found)
    # print('Predicted: ', num_tool)
    # print('Found: ', num_true)

    precision = float(100 * num_true / num_tool)
    recall = float(100 * num_true / num_ref)
    fScore = float(2 * precision*recall/(precision+recall))
    print('P: %.3f R: %.3f F: %.3f' % (precision, recall, fScore))

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

    return dis, found, sup


def bpePlot(xlabel, bpe):
    fig, ax = plt.subplots()
    # ax.set_xlabel(xlabel)
    ax.set_ylabel('Breakpoint error (bp)')
    ax.boxplot(bpe, showfliers=False)
    ax.set_xticklabels(['SoftSV', 'Lumpy', 'TIDDIT', 'OUR'])


def supPlot(sup):
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    bins_sz = [2*x for x in range(25)]

    ax1.hist(sup[0], bins=bins_sz, ec='black')
    ax1.set_ylim([0, 500])
    ax1.set_ylabel('Number of deletions')
    ax1.set_xlabel('Breakpoint Support (SC)')
    ax1.set_title('SoftSV')
    ax1.axvline(median(sup[0]), linestyle='dashed',
                color='k', linewidth=1)

    ax2.hist(sup[1], bins=bins_sz, ec='black')
    ax2.set_ylim([0, 500])
    ax2.set_ylabel('Number of deletions')
    ax2.set_xlabel('Breakpoint Support (SR)')
    ax2.set_title('Lumpy')
    ax2.axvline(median(sup[1]), linestyle='dashed',
                color='k', linewidth=1)

    ax3.hist(sup[2], bins=bins_sz, ec='black')
    ax3.set_ylim([0, 500])
    ax3.set_ylabel('Number of deletions')
    ax3.set_xlabel('Breakpoint Support (SR)')
    ax3.set_title('TIDDIT')
    ax3.axvline(median(sup[2]), linestyle='dashed',
                color='k', linewidth=1)

    ax4.hist(sup[3], bins=bins_sz, ec='black')
    ax4.set_ylim([0, 500])
    ax4.set_ylabel('Number of deletions')
    ax4.set_xlabel('Breakpoint Support (SC + SR)')
    ax4.set_title('OUR')
    ax4.axvline(median(sup[3]), linestyle='dashed',
                color='k', linewidth=1)

    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    # fig.savefig('filename.eps', format='eps')


def real():
    # Benchmark
    real_dels_svclassify_ref = readSvclassify(
        '/Users/alok/IIIT/GS/Personalis_1000_Genomes_deduplicated_deletions.bed', True)
    # real_ins_svclassify_ref = readSvclassify(
    #     '/Users/alok/IIIT/GS/Spiral_Genetics_insertions.bed', False)

    # Tools
    # ----------------------------------------------------------------------
    # Lumpy
    # -----------------------------------
    real_dels_lumpy = readLumpy(
        '/Users/alok/IIIT/Real/Results/lumpy/lumpy_real.vcf')
    bpe_lumpy_real, f_lumpy_real, sup_lumpy_real = compare(real_dels_lumpy, real_dels_svclassify_ref,
                                                           'REAL-DELS-LUMPY', True, RO=REAL_RO)
    # ----------------------------------------------------------------------
    # Tiddit
    # -----------------------------------
    real_dels_tiddit = readTiddit(
        '/Users/alok/IIIT/Real/Results/tiddit/tiddit_real.vcf')
    bpe_tiddit_real, f_tiddit_real, sup_tiddit_real = compare(real_dels_tiddit, real_dels_svclassify_ref,
                                                              'REAL-DELS-TIDDIT', True, RO=REAL_RO)
    # ----------------------------------------------------------------------
    # Softsv
    # -----------------------------------
    real_dels_softsv = readSoftsv(
        '/Users/alok/IIIT/Real/Results/softsv/deletions_small.txt', '/Users/alok/IIIT/Real/Results/softsv/deletions.txt')
    bpe_softsv_real, f_softsv_real, sup_softsv_real = compare(real_dels_softsv, real_dels_svclassify_ref,
                                                              'REAL-DELS-SOFTSV', True, RO=REAL_RO)
    # ----------------------------------------------------------------------
    # Minus
    # -----------------------------------
    real_dels_minus = readMinus(
        '/Users/alok/IIIT/Real/Results/plusminus/output.vcf')
    bpe_minus_real, f_minus_real, sup_minus_real = compare(real_dels_minus, real_dels_svclassify_ref,
                                                           'REAL-DELS-MINUS', True, RO=REAL_RO)

    # Plot breakpoint support
    if SUPPORT_PLOT:
        sup_real = [sup_softsv_real, sup_lumpy_real,
                    sup_tiddit_real, sup_minus_real]
        supPlot(sup_real)
        plt.show()


def sim():
    # Benchmark
    sim_dels_ref = readSim('/Users/alok/IIIT/Simulations/sim1.txt', True)

    # Tools
    # ----------------------------------------------------------------------
    # Lumpy
    # -----------------------------------
    sim_dels_lumpy_5x = readLumpy(
        '/Users/alok/IIIT/Simulations/Results/5x/lumpy/lumpy_sim_5x.vcf')
    sim_dels_lumpy_10x = readLumpy(
        '/Users/alok/IIIT/Simulations/Results/10x/lumpy/lumpy_sim_10x.vcf')
    sim_dels_lumpy_20x = readLumpy(
        '/Users/alok/IIIT/Simulations/Results/20x/lumpy/lumpy_sim_20x.vcf')
    sim_dels_lumpy_30x = readLumpy(
        '/Users/alok/IIIT/Simulations/Results/30x/lumpy/lumpy_sim_30x.vcf')
    # -----------------------------------
    compare(sim_dels_lumpy_5x, sim_dels_ref,
            'SIM-DELS-LUMPY-5x', True, RO=SIM_RO)
    compare(sim_dels_lumpy_10x, sim_dels_ref,
            'SIM-DELS-LUMPY-10x', True, RO=SIM_RO)
    compare(sim_dels_lumpy_20x, sim_dels_ref,
            'SIM-DELS-LUMPY-20x', True, RO=SIM_RO)
    bpe_lumpy_30x, f_lumpy_30x, sup_lumpy_30x = compare(
        sim_dels_lumpy_30x, sim_dels_ref, 'SIM-DELS-LUMPY-30x', True, RO=SIM_RO)
    # ----------------------------------------------------------------------
    # Tiddit
    # -----------------------------------
    sim_dels_tiddit_5x = readTiddit(
        '/Users/alok/IIIT/Simulations/Results/5x/tiddit/tiddit_sim_5x.vcf')
    sim_dels_tiddit_10x = readTiddit(
        '/Users/alok/IIIT/Simulations/Results/10x/tiddit/tiddit_sim_10x.vcf')
    sim_dels_tiddit_20x = readTiddit(
        '/Users/alok/IIIT/Simulations/Results/20x/tiddit/tiddit_sim_20x.vcf')
    sim_dels_tiddit_30x = readTiddit(
        '/Users/alok/IIIT/Simulations/Results/30x/tiddit/tiddit_sim_30x.vcf')
    # -----------------------------------
    compare(sim_dels_tiddit_5x, sim_dels_ref,
            'SIM-DELS-TIDDIT-5x', True, RO=SIM_RO)
    compare(sim_dels_tiddit_10x, sim_dels_ref,
            'SIM-DELS-TIDDIT-10x', True, RO=SIM_RO)
    compare(sim_dels_tiddit_20x, sim_dels_ref,
            'SIM-DELS-TIDDIT-20x', True, RO=SIM_RO)
    bpe_tiddit_30x, f_tiddit_30x, sup_tiddit_30x = compare(
        sim_dels_tiddit_30x, sim_dels_ref, 'SIM-DELS-TIDDIT-30x', True, RO=SIM_RO)
    # ----------------------------------------------------------------------
    # Softsv
    # -----------------------------------
    sim_dels_softsv_5x = readSoftsv('/Users/alok/IIIT/Simulations/Results/5x/softsv/deletions_small.txt',
                                    '/Users/alok/IIIT/Simulations/Results/5x/softsv/deletions.txt')
    sim_dels_softsv_10x = readSoftsv('/Users/alok/IIIT/Simulations/Results/10x/softsv/deletions_small.txt',
                                     '/Users/alok/IIIT/Simulations/Results/10x/softsv/deletions.txt')
    sim_dels_softsv_20x = readSoftsv('/Users/alok/IIIT/Simulations/Results/20x/softsv/deletions_small.txt',
                                     '/Users/alok/IIIT/Simulations/Results/20x/softsv/deletions.txt')
    sim_dels_softsv_30x = readSoftsv('/Users/alok/IIIT/Simulations/Results/30x/softsv/deletions_small.txt',
                                     '/Users/alok/IIIT/Simulations/Results/30x/softsv/deletions.txt')
    # -----------------------------------
    compare(sim_dels_softsv_5x, sim_dels_ref,
            'SIM-DELS-SOFTSV-5x', True, RO=SIM_RO)
    compare(sim_dels_softsv_10x, sim_dels_ref,
            'SIM-DELS-SOFTSV-10x', True, RO=SIM_RO)
    compare(sim_dels_softsv_20x, sim_dels_ref,
            'SIM-DELS-SOFTSV-20x', True, RO=SIM_RO)
    bpe_softsv_30x, f_softsv_30x, sup_softsv_30x = compare(
        sim_dels_softsv_30x, sim_dels_ref, 'SIM-DELS-SOFTSV-30x', True, RO=SIM_RO)

    # ----------------------------------------------------------------------
    # Minus
    # -----------------------------------
    sim_dels_minus_5x = readMinus(
        '/Users/alok/IIIT/Simulations/Results/5x/plusminus/output.vcf')
    sim_dels_minus_10x = readMinus(
        '/Users/alok/IIIT/Simulations/Results/10x/plusminus/output.vcf')
    sim_dels_minus_20x = readMinus(
        '/Users/alok/IIIT/Simulations/Results/20x/plusminus/output.vcf')
    sim_dels_minus_30x = readMinus(
        '/Users/alok/IIIT/Simulations/Results/30x/plusminus/output.vcf')
    # -----------------------------------
    compare(sim_dels_minus_5x, sim_dels_ref,
            'SIM-DELS-MINUS-5x', True, RO=SIM_RO)
    compare(sim_dels_minus_10x, sim_dels_ref,
            'SIM-DELS-MINUS-10x', True, RO=SIM_RO)
    compare(sim_dels_minus_20x, sim_dels_ref,
            'SIM-DELS-MINUS-20x', True, RO=SIM_RO)
    bpe_minus_30x, f_minus_30x, sup_minus_30x = compare(
        sim_dels_minus_30x, sim_dels_ref, 'SIM-DELS-MINUS-30x', True, RO=SIM_RO)

    # Plot breakpoint error
    if BREAKPOINT_PLOT:
        bpe_30x = [bpe_softsv_30x, bpe_lumpy_30x,
                   bpe_tiddit_30x, bpe_minus_30x]
        bpePlot('All', bpe_30x)
        plt.show()

    # Plot breakpoint support
    if SUPPORT_PLOT:
        sup_30x = [sup_softsv_30x, sup_lumpy_30x,
                   sup_tiddit_30x, sup_minus_30x]
        supPlot(sup_30x)
        plt.show()


def main():
    real()
    sim()
    return


if __name__ == '__main__':
    main()
