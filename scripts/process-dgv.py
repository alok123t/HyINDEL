with open('/Users/alok/Data/NA12878/GRCh37_hg19_variants_2014-10-16.txt') as inp, open('/Users/alok/Data/NA12878/dgv-full.bed', 'w') as out:
    chrz = {}
    var = {}
    sub = {}
    seq = {}
    for line in inp:
        if line.find('NA12878') == -1:
            continue
        l = line.split()
        if l[0] == 'variantaccession':
            continue
        if l[5] != 'deletion':
            continue
        sz = int(l[3])- int(l[2]) + 1
        if not(sz >= 200 and sz <= 10000):
            continue
        out.write(line)
        # out.write('chr' + l[1] + '\t' + l[2] + '\t' + l[3] + '\t' + str(sz) + '\n')

        # if l[1] not in chrz:
        #     chrz[l[1]] = 0
        # if l[4] not in var:
        #     var[l[4]] = 0
        # if l[5] not in sub:
        #     sub[l[5]] = 0
        # if l[8] not in seq:
        #     seq[l[8]] = 0
        # chrz[l[1]] += 1
        # var[l[4]] += 1
        # sub[l[5]] += 1
        # seq[l[8]] += 1
    # print(chrz)
    # print(var)
    # print(sub)
    # print(seq)

# variantaccession	chr	start	end	varianttype	variantsubtype	reference	pubmedid	method	platform	mergedvariants	supportingvariants	mergedorsample	frequency	samplesize	observedgains	observedlosses	cohortdescription	genes	samples
# nsv482937	1	1	2300000	CNV	loss	Iafrate_et_al_2004	15286789	BAC aCGH,FISH			nssv2995976	M		39	0	1		ACAP3,AGRN,ANKRD65,ATAD3A,ATAD3B,ATAD3C,AURKAIP1,B3GALT6,C1orf159,C1orf170,C1orf233,C1orf86,CALML6,CCNL2,CDK11A,CDK11B,CPSF3L,DDX11L1,DVL1,FAM132A,FAM138A,FAM138F,FAM41C,FAM87B,GABRD,GLTPD1,GNB1,HES4,ISG15,KIAA1751,KLHL17,LINC00115,LINC01128,LOC100129534,LOC100130417,LOC100132062,LOC100132287,LOC100133331,LOC100288069,LOC148413,LOC254099,LOC729737,MIB2,MIR200A,MIR200B,MIR429,MIR6723,MIR6726,MIR6727,MIR6808,MIR6859-1,MIR6859-2,MMP23A,MMP23B,MORN1,MRPL20,MXRA8,NADK,NOC2L,OR4F16,OR4F29,OR4F3,OR4F5,PLEKHN1,PRKCZ,PUSL1,RNF223,SAMD11,SCNN1D,SDF4,SKI,SLC35E2,SLC35E2B,SSU72,TAS1R3,TMEM240,TMEM52,TMEM88B,TNFRSF18,TNFRSF4,TTLL10,UBE2J2,VWA1,WASH7P
