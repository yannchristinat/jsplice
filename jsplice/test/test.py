from jsplice.lib.bio import Region
from jsplice.lib.asm import ASM
from jsplice.lib.expdesign import AbstractDesign

j1 = Region('chr13', 107029433, 107029547, '+')
j2 = Region('chr13', 107029433, 107030255, '+')
j3 = Region('chr13', 107029433, 107030363, '+')

jxnCnts = dict()
jxnCnts[j1] = AbstractDesign()
jxnCnts[j1].addElement('rep2','KD', 148)
jxnCnts[j1].addElement('rep2','ctrl', 134)
jxnCnts[j1].addElement('rep1','KD', 83)
jxnCnts[j1].addElement('rep1','ctrl', 104)

jxnCnts[j2] = AbstractDesign()
jxnCnts[j2].addElement('rep2','KD', 8)
jxnCnts[j2].addElement('rep2','ctrl', 8)
jxnCnts[j2].addElement('rep1','KD', 15)
jxnCnts[j2].addElement('rep1','ctrl', 15)

jxnCnts[j3] = AbstractDesign()
jxnCnts[j3].addElement('rep2','KD', 7)
jxnCnts[j3].addElement('rep2','ctrl', 8)
jxnCnts[j3].addElement('rep1','KD', 11)
jxnCnts[j3].addElement('rep1','ctrl', 9)

expdesign = AbstractDesign()
expdesign.addElement('rep1', 'KD', 0)
expdesign.addElement('rep1', 'ctrl', 0)
expdesign.addElement('rep2', 'KD', 0)
expdesign.addElement('rep2', 'ctrl', 0)
expdesign.cond0 = 'KD'
expdesign.cond1 = 'ctrl'

a = ASM([j1, j2, j3])
a.exons = []
a.assignCounts(jxnCnts, dict(), expdesign)
rdiffJxn, pvalJxn = a.getLargestRDiff(ASM.JXN, 2, 0.1, 2, 0.05, 1.5)

print rdiffJxn
print pvalJxn

