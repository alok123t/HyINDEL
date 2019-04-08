from statistics import median
import numpy as np
import matplotlib.pyplot as plt

# Remove events < 50 bp in tools
FLAG_50 = True
# Plot breakpoint error plot
BREAKPOINT_PLOT = False
# Plot support
SUPPORT_PLOT = False
# TODO: Implement seperate results for small and large variants
# Print metrics seperately for large and small variants
PRINT_SPLIT = False
SPLIT_LARGE = 500
# Maximum error in breakpoint for insertions
DIST_INS = 10
# TODO: Implement sequence comparison
# Compare insertion sequence
COMPARE_INS_SEQ = False
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

    if COMPARE_INS_SEQ:
        f_seq = False
        if abs(len(a.seq) - len(b.seq)) <= DIST_INS:
            f_seq = True
        return f_pos & f_seq
    else:
        return f_pos


class Deletion:
    def __init__(self, t_chr, t_st, t_en, t_iD=-1, t_pe=-1, t_sc=-1, t_sr=-1):
        self.chr = modChr(t_chr)
        self.st = t_st
        self.en = t_en
        self.iD = t_iD
        self.pe = t_pe
        self.sc = t_sc
        self.sr = t_sr

    def __str__(self):
        return self.chr+'\t' + str(self.st)+'\t' + str(self.en)


class Insertion:
    def __init__(self, t_chr, t_pos, t_seq='', t_seqLen=-1, t_refPos=-1, t_iD=-1, t_sc=-1):
        self.chr = modChr(t_chr)
        self.pos = t_pos  # Position on reference
        self.seq = t_seq
        self.seqLen = t_seqLen
        self.refPos = t_refPos  # Position on sample
        self.iD = t_iD
        self.sc = t_sc

    def __str__(self):
        return self.chr+'\t' + str(self.pos)


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
            ret.append(Deletion(t_chr=l[0], t_st=int(l[1]), t_en=int(l[2])))
        else:
            ret.append(Insertion(t_chr=l[0], t_pos=int(l[1])))
            ret.append(Insertion(t_chr=l[0], t_pos=int(l[2])))

    return ret


def readMy2(fName):
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
        sup_pe = int(format_col[1])
        sup_sr = int(format_col[2])
        sup_sc = int(format_col[3])

        ret.append(Deletion(t_chr=l[0], t_st=int(
            l[1]), t_en=int(info_col[2].split('=')[1]), t_pe=sup_pe, t_sc=sup_sc, t_sr=sup_sr))

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
        if l[6] != 'PASS':
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


def readSim(bedpeFName, eventFName, fastaFName, delsFlag):
    inp = open(bedpeFName)
    lines = inp.readlines()
    inp.close()

    ret = []
    for line in lines:
        l = line.split()
        if l[0] == 'LITERAL':
            continue
        if delsFlag:
            if 'DEL' in l[6]:
                ret.append(
                    Deletion(t_chr=l[0], t_st=int(l[2]), t_en=int(l[4])))
        else:
            if 'INC' in l[6]:
                ret.append(Insertion(t_chr=l[0], t_pos=int(
                    l[2]), t_iD=l[6].split('::')[1]))

    if not delsFlag:
        inp = open(eventFName)
        lines = inp.readlines()
        inp.close()

        for line in lines:
            l = line.split()
            if l[1] == 'INC':
                t_iD = l[0]
                for i in range(len(ret)):
                    if ret[i].iD == modChr(t_iD):
                        ret[i].seqLen = int(l[2])
                        ret[i].refPos = int(l[4])
                        break

    return ret


"""
returns [[breakpoint error], [found], [support]]
[breakpoint error]: only for detected variants
[found]: True/False if present in benchmark
[support]: sc/sr support only for detected variants (if present)
"""


def compare(toolUnfiltered, ref, toolName, delsFlag):
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
    for i in range(len(tool)):
        ok = False
        for j in range(len(ref)):
            if delsFlag:
                if checkDel(tool[i], ref[j]):
                    # if found[j]:
                    #     print('Double:', tool[i])
                    found[j] = True
                    ok = True
                    dis.append(abs(tool[i].en - ref[j].en) +
                               abs(tool[i].st - ref[j].st))
                    exact_sup = tool[i].sc + tool[i].sr
                    if exact_sup > 0:
                        sup.append(exact_sup)
            else:
                if checkIns(tool[i], ref[j]):
                    found[j] = True
                    dis.append(abs(tool[i].pos - ref[j].pos))
                    print(dis[-1])
        # if not ok:
        #     print('Miss', tool[i])

    num_ref = len(ref)
    num_tool = len(tool)
    num_true = sum(found)
    # print('Predicted: ', num_tool)
    # print('Found: ', num_true)

    precision = float(100 * num_true / num_tool)
    recall = float(100*num_true / num_ref)
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
        print('Large P: %.3f R: %.3f F: %.3f' %
              (split_pr_large, split_re_large, split_f_large))
        split_pr_small = float(100 * num_small / tool_small)
        split_re_small = float(100 * num_small / ref_small)
        split_f_small = float(
            2 * split_pr_small * split_re_small / (split_pr_small + split_re_small))
        print('Small P: %.3f R: %.3f F: %.3f' %
              (split_pr_small, split_re_small, split_f_small))
    print('-' * 27)

    return dis, found, sup


def find_missing(a, b, ref):
    pre_a = []
    pre_b = []
    print('A: ', sum(a))
    print('B: ', sum(b))
    for i in range(len(a)):
        if a[i] & (not b[i]):
            pre_a.append(i)
        elif (not a[i]) & b[i]:
            pre_b.append(i)
    print('Only in A: ', len(pre_a))
    print('Only in B: ', len(pre_b))
    for x in pre_b:
        if x not in pre_a:
            print(ref[x].chr, ref[x].st, ref[x].en, ref[x].en-ref[x].st+1)


def real():
    # Benchmark
    real_dels_svclassify_ref = readSvclassify(
        '/Users/alok/Data/GS/Personalis_1000_Genomes_deduplicated_deletions.bed', True)
    real_ins_svclassify_ref = readSvclassify(
        '/Users/alok/Data/GS/Spiral_Genetics_insertions.bed', False)

    # Tools
    # ----------------------------------------------------------------------
    # Lumpy
    # -----------------------------------
    # real_dels_lumpy = readLumpy(
    #     '/Users/alok/Data/Results/Real/lumpy/lumpy_real.vcf')
    # compare(real_dels_lumpy, real_dels_svclassify_ref, 'REAL-DELS-LUMPY', True)
    # ----------------------------------------------------------------------
    # Tiddit
    # -----------------------------------
    # real_dels_tiddit = readTiddit(
    #     '/Users/alok/Data/Results/Real/tiddit/tiddit_real.vcf')
    # compare(real_dels_tiddit, real_dels_svclassify_ref,
    # 'REAL-DELS-TIDDIT', True)
    # ----------------------------------------------------------------------
    # Softsv
    # -----------------------------------
    # real_dels_softsv = readSoftsv(
    #     '/Users/alok/Data/Results/Real/softsv/deletions_small.txt', '/Users/alok/Data/Results/Real/softsv/deletions.txt')
    # compare(real_dels_softsv, real_dels_svclassify_ref,
    # 'REAL-DELS-SOFTSV', True)
    # ----------------------------------------------------------------------
    # My deletions
    # -----------------------------------
    # real_my_ins = readMy('/Users/alok/tmp/today/my/tmp/ins.txt', False)
    # compare(real_my_ins, real_ins_svclassify_ref, 'REAL-INS-MY', False)
    # return
    real_my_latest = readMy2('/Users/alok/tmp/April/7/my/output.vcf')
    bpe_my_latest, f_my_latest, sup_my_latest = compare(real_my_latest, real_dels_svclassify_ref,
                                                        'REAL-DELS-MY-LATEST', True)
    real_dels_my_flanks_3 = readMy2(
        '/Users/alok/tmp/old/2019/Mar/29/my3/output.vcf')
    bpe_my_prv, f_my_prv, sup_my_prv = compare(real_dels_my_flanks_3, real_dels_svclassify_ref,
                                               'REAL-DELS-MY-FLANKS-10,20', True)
    # find_missing(f_my_latest, f_my_prv, real_dels_svclassify_ref)
    return
    real_dels_my_posttry_3 = readMy2('/Users/alok/tmp/today/output.vcf')
    compare(real_dels_my_posttry_3, real_dels_svclassify_ref,
            'REAL-DELS-MY-POSTTRY-10,20', True)
    real_dels_my_posttry_2 = readMy2(
        '/Users/alok/Data/Results/Real/my_posttry/output.vcf')
    compare(real_dels_my_posttry_2, real_dels_svclassify_ref,
            'REAL-DELS-MY-POSTTRY-15,30', True)
    real_dels_my_flanks = readMy2('/Users/alok/tmp/today/my/output.vcf')
    compare(real_dels_my_flanks, real_dels_svclassify_ref,
            'REAL-DELS-MY-FLANKS', True)
    real_dels_my_flanks_2 = readMy2('/Users/alok/tmp/today/my2/output.vcf')
    compare(real_dels_my_flanks_2, real_dels_svclassify_ref,
            'REAL-DELS-MY-FLANKS-15,30', True)
    # real_dels_fixpost = readMy2('/Users/alok/tmp/today/flank/output.vcf')
    # compare(real_dels_fixpost, real_dels_svclassify_ref,
    #         'REAL-DELS-MY-FIXPOST', True)
    # real_dels_my_preds_all = readMy('/Users/alok/tmp/today/preds_all.bed', True)
    # compare(real_dels_my_preds_all, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-ALL', True)
    # real_dels_my_preds_flanks = readMy('/Users/alok/tmp/today/preds_flanks.bed', True)
    # compare(real_dels_my_preds_flanks, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-FLANKS', True)
    # real_dels_my_preds_sup_10 = readMy('/Users/alok/tmp/today/preds_sup_10.bed', True)
    # compare(real_dels_my_preds_sup_10, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-10', True)
    # real_dels_my_preds_sup_15 = readMy('/Users/alok/tmp/today/preds_sup_15.bed', True)
    # compare(real_dels_my_preds_sup_15, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-15', True)
    # real_dels_my_preds_sup_10_20 = readMy('/Users/alok/tmp/today/preds_sup_10_20.bed', True)
    # compare(real_dels_my_preds_sup_10_20, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-SEPERATE-10-20', True)
    # real_dels_my_preds_sup_10_20_flanks = readMy('/Users/alok/tmp/today/preds_sup_10_20_flanks.bed', True)
    # compare(real_dels_my_preds_sup_10_20_flanks, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-SEPERATE-10-20-FLANKS', True)
    # real_dels_my_preds_sup_15_30 = readMy('/Users/alok/tmp/today/preds_sup_15_30.bed', True)
    # compare(real_dels_my_preds_sup_15_30, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-SEPERATE-15-30', True)
    # real_dels_my_preds_sup_15_30_flanks = readMy('/Users/alok/tmp/today/preds_sup_15_30_flanks.bed', True)
    # compare(real_dels_my_preds_sup_15_30_flanks, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-SEPERATE-15-30-FLANKS', True)
    # real_dels_my_preds_sup_15_30_flanks_5000 = readMy('/Users/alok/tmp/today/preds_sup_15_30_flanks_5000.bed', True)
    # compare(real_dels_my_preds_sup_15_30_flanks_5000, real_dels_svclassify_ref,
    #     'REAL-DELS-MY-PREDS-SUP-SEPERATE-15-30-FLANKS-5000', True)
    # ----------------------------------------------------------------------

    # real_dels_cov_fail = readMy('/Users/alok/tmp/today/cov_fail.bed', True)
    # real_dels_fixpost_bed = readMy('/Users/alok/tmp/today/output_all.bed', True)
    # real_dels_fixpost_merge_bed = readMy('/Users/alok/tmp/today/merge_all.bed', True)
    # real_dels_my_prv = readMy(
    #     '/Users/alok/Data/Results/Real/my/deletions.bed', True)
    # real_dels_my = readMy2('/Users/alok/Data/Results/Real/my/output.vcf')
    # real_dels_my_1 = readMy('/Users/alok/Tools/indel-detect/build/1.bed', True)
    # real_dels_my_1_merge = readMy(
    #     '/Users/alok/Tools/indel-detect/build/1_merge.bed', True)
    # real_dels_my_2 = readMy('/Users/alok/Tools/indel-detect/build/2.bed', True)
    # real_dels_my_2_merge = readMy(
    #     '/Users/alok/Tools/indel-detect/build/2_merge.bed', True)

    # dis_my_prv, found_my_prv = compare(real_dels_my_prv, real_dels_svclassify_ref,
    #                                    'REAL-DELS-MY-PRV', True)
    # compare(real_dels_fixpost_bed, real_dels_svclassify_ref,
    #         'REAL-DELS-MY-FIXPOST-BED', True)
    # compare(real_dels_fixpost_merge_bed, real_dels_svclassify_ref,
    #         'REAL-DELS-MY-FIXPOST-MERGE-BED', True)
    # compare(real_dels_my, real_dels_svclassify_ref,
    #    'REAL-DELS-MY-RAW', True)
    # dis_my, found_my = compare(real_dels_my_1, real_dels_svclassify_ref, 'REAL-DELS-MY-1', True)
    # dis_my, found_my = compare(real_dels_my_1_merge, real_dels_svclassify_ref, 'REAL-DELS-MY-1-MERGE', True)
    # dis_my, found_my = compare(real_dels_my_2, real_dels_svclassify_ref, 'REAL-DELS-MY-2', True)
    # dis_my, found_my= compare(real_dels_my_2_merge, real_dels_svclassify_ref,'REAL-DELS-MY-2-MERGE', True)
    # compare(real_dels_manta, real_dels_svclassify_ref,'REAL-DELS-MANTA', True)
    # find_missing(found_lumpy, found_my, real_dels_svclassify_ref)
    # find_missing(found_tiddit, found_my)
    # find_missing(found_manta, found_my, real_dels_svclassify_ref)

    # Plot breakpoint error
    # if BREAKPOINT_PLOT:
    #     dis=[dis_lumpy, dis_tiddit, dis_softsv, dis_my]
    #     fig, ax=plt.subplots()
    #     ax.set_ylabel('Breakpoint error (bp)')
    #     ax.boxplot(dis, showfliers = False)
    #     plt.show()


def sim():
    # Benchmark
    sim_dels_ref = readSim('/Users/alok/Data/Simulations/sim_ref.bedpe',
                           '/Users/alok/Data/Simulations/sim_ref.event', '/Users/alok/Data/Simulations/sim_ref.fasta', True)
    # sim_ins_ref = readSim('/Users/alok/Data/Simulations/sim_ref.bedpe',
    #                       '/Users/alok/Data/Simulations/sim_ref.event', '/Users/alok/Data/Simulations/sim_ref.fasta', False)
    """
    # Tools
    # ----------------------------------------------------------------------
    # Lumpy
    # -----------------------------------
    sim_dels_lumpy_5x = readLumpy(
        '/Users/alok/Data/Results/Simulations/5x/lumpy/lumpy_sim_5x.vcf')
    sim_dels_lumpy_10x = readLumpy(
        '/Users/alok/Data/Results/Simulations/10x/lumpy/lumpy_sim_10x.vcf')
    sim_dels_lumpy_20x = readLumpy(
        '/Users/alok/Data/Results/Simulations/20x/lumpy/lumpy_sim_20x.vcf')
    sim_dels_lumpy_30x = readLumpy(
        '/Users/alok/Data/Results/Simulations/30x/lumpy/lumpy_sim_30x.vcf')
    # -----------------------------------
    compare(sim_dels_lumpy_5x, sim_dels_ref, 'SIM-DELS-LUMPY-5x', True)
    compare(sim_dels_lumpy_10x, sim_dels_ref, 'SIM-DELS-LUMPY-10x', True)
    compare(sim_dels_lumpy_20x, sim_dels_ref, 'SIM-DELS-LUMPY-20x', True)
    bpe_lumpy_30x, f_lumpy_30x, sup_lumpy_30x = compare(
        sim_dels_lumpy_30x, sim_dels_ref, 'SIM-DELS-LUMPY-30x', True)
    # ----------------------------------------------------------------------
    # Tiddit
    # -----------------------------------
    sim_dels_tiddit_5x = readTiddit(
        '/Users/alok/Data/Results/Simulations/5x/tiddit/tiddit_sim_5x.vcf')
    sim_dels_tiddit_10x = readTiddit(
        '/Users/alok/Data/Results/Simulations/10x/tiddit/tiddit_sim_10x.vcf')
    sim_dels_tiddit_20x = readTiddit(
        '/Users/alok/Data/Results/Simulations/20x/tiddit/tiddit_sim_20x.vcf')
    sim_dels_tiddit_30x = readTiddit(
        '/Users/alok/Data/Results/Simulations/30x/tiddit/tiddit_sim_30x.vcf')
    # -----------------------------------
    compare(sim_dels_tiddit_5x, sim_dels_ref, 'SIM-DELS-TIDDIT-5x', True)
    compare(sim_dels_tiddit_10x, sim_dels_ref, 'SIM-DELS-TIDDIT-10x', True)
    compare(sim_dels_tiddit_20x, sim_dels_ref, 'SIM-DELS-TIDDIT-20x', True)
    bpe_tiddit_30x, f_tiddit_30x, sup_tiddit_30x = compare(
        sim_dels_tiddit_30x, sim_dels_ref, 'SIM-DELS-TIDDIT-30x', True)
    # ----------------------------------------------------------------------
    # Softsv
    # -----------------------------------
    sim_dels_softsv_5x = readSoftsv('/Users/alok/Data/Results/Simulations/5x/softsv/deletions_small.txt',
                                    '/Users/alok/Data/Results/Simulations/5x/softsv/deletions.txt')
    sim_dels_softsv_10x = readSoftsv('/Users/alok/Data/Results/Simulations/10x/softsv/deletions_small.txt',
                                     '/Users/alok/Data/Results/Simulations/10x/softsv/deletions.txt')
    sim_dels_softsv_20x = readSoftsv('/Users/alok/Data/Results/Simulations/20x/softsv/deletions_small.txt',
                                     '/Users/alok/Data/Results/Simulations/20x/softsv/deletions.txt')
    sim_dels_softsv_30x = readSoftsv('/Users/alok/Data/Results/Simulations/30x/softsv/deletions_small.txt',
                                     '/Users/alok/Data/Results/Simulations/30x/softsv/deletions.txt')
    # -----------------------------------
    compare(sim_dels_softsv_5x, sim_dels_ref, 'SIM-DELS-SOFTSV-5x', True)
    compare(sim_dels_softsv_10x, sim_dels_ref, 'SIM-DELS-SOFTSV-10x', True)
    compare(sim_dels_softsv_20x, sim_dels_ref, 'SIM-DELS-SOFTSV-20x', True)
    bpe_softsv_30x, f_softsv_30x, sup_softsv_30x = compare(
        sim_dels_softsv_30x, sim_dels_ref, 'SIM-DELS-SOFTSV-30x', True)
    """
    # ----------------------------------------------------------------------
    # My deletions
    # -----------------------------------
    sim_dels_my_5x = readMy2(
        '/Users/alok/Data/Results/Simulations/5x/my/output.vcf')
    sim_dels_my_10x = readMy2(
        '/Users/alok/Data/Results/Simulations/10x/my/output.vcf')
    sim_dels_my_20x = readMy2(
        '/Users/alok/Data/Results/Simulations/20x/my/output.vcf')
    sim_dels_my_30x = readMy2(
        '/Users/alok/Data/Results/Simulations/20x/my/output.vcf')
    # -----------------------------------
    compare(sim_dels_my_5x, sim_dels_ref, 'SIM-DELS-MY-5x', True)
    compare(sim_dels_my_10x, sim_dels_ref, 'SIM-DELS-MY-10x', True)
    compare(sim_dels_my_20x, sim_dels_ref, 'SIM-DELS-MY-20x', True)
    bpe_my_30x, f_my_30x, sup_my_30x = compare(
        sim_dels_my_30x, sim_dels_ref, 'SIM-DELS-MY-30x', True)
    """
    # ----------------------------------------------------------------------
    # My insertions
    # -----------------------------------
    sim_ins_my_10x = readMy(
        '/Users/alok/tmp/today/my/10x/insertions.bed', False)
    sim_ins_my_20x = readMy(
        '/Users/alok/tmp/today/my/20x/insertions.bed', False)
    sim_ins_my_30x = readMy(
        '/Users/alok/tmp/today/my/30x/insertions.bed', False)
    # -----------------------------------
    compare(sim_ins_my_10x, sim_ins_ref, 'SIM-INS-MY-10x', False)
    compare(sim_ins_my_20x, sim_ins_ref, 'SIM-INS-MY-20x', False)
    # compare_ins(sim_ins_my_30x, sim_ins_ref, 'SIM-INS-MY-30x', False)
    # ----------------------------------------------------------------------
    """

    # Plot breakpoint error
    if BREAKPOINT_PLOT:

        def simBpePlot(xlabel, bpe):
            fig, ax = plt.subplots()
            # ax.set_xlabel(xlabel)
            ax.set_ylabel('Breakpoint error (bp)')
            ax.boxplot(bpe, showfliers=False)
            ax.set_xticklabels(['SoftSV', 'Lumpy', 'TIDDIT', 'OUR'])

        bpe_30x = [bpe_softsv_30x, bpe_lumpy_30x, bpe_tiddit_30x, bpe_my_30x]

        simBpePlot('All', bpe_30x)

        plt.show()

    # Plot supports
    if SUPPORT_PLOT:

        def simSupPlot(sup):
            fig = plt.figure()

            ax1 = fig.add_subplot(2, 2, 1)
            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4)

            bins_sz = [2*x for x in range(25)]

            ax1.hist(sup[0], bins=bins_sz, ec='black')
            ax1.set_ylim([0, 500])
            ax1.set_ylabel('Number of deletions')
            ax1.set_xlabel('Support (SC)')
            ax1.set_title('SoftSV')
            ax1.axvline(median(sup[0]), linestyle='dashed',
                        color='k', linewidth=1)

            ax2.hist(sup[1], bins=bins_sz, ec='black')
            ax2.set_ylim([0, 500])
            ax2.set_ylabel('Number of deletions')
            ax2.set_xlabel('Support (SR)')
            ax2.set_title('Lumpy')
            ax2.axvline(median(sup[1]), linestyle='dashed',
                        color='k', linewidth=1)

            ax3.hist(sup[2], bins=bins_sz, ec='black')
            ax3.set_ylim([0, 500])
            ax3.set_ylabel('Number of deletions')
            ax3.set_xlabel('Support (SR)')
            ax3.set_title('TIDDIT')
            ax3.axvline(median(sup[2]), linestyle='dashed',
                        color='k', linewidth=1)

            ax4.hist(sup[3], bins=bins_sz, ec='black')
            ax4.set_ylim([0, 500])
            ax4.set_ylabel('Number of deletions')
            ax4.set_xlabel('Support (SC + SR)')
            ax4.set_title('OUR')
            ax4.axvline(median(sup[3]), linestyle='dashed',
                        color='k', linewidth=1)

            plt.subplots_adjust(wspace=0.5, hspace=0.5)
            # fig.savefig('filename.eps', format='eps')

        sup_30x = [sup_softsv_30x, sup_lumpy_30x, sup_tiddit_30x, sup_my_30x]

        simSupPlot(sup_30x)

        plt.show()


def main():
    real()
    # sim()
    return


if __name__ == '__main__':
    main()
