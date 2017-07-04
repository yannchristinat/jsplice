'''This file is part of jSplice.

    jSplice is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jSplice is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jSplice.  If not, see <http://www.gnu.org/licenses/>
    
    (C) copyright ETH Zurich; 2014 Institute for Molecular Health Sciences, Krek group
    Developer: Yann Christinat
    
'''

import bio, math
from scipy.stats import fisher_exact

# Module describing alternatively spliced modules or ASMs
class ASM:
    JXN=1
    EXN=2
    
    def __init__(self, junctions):
        self.junctions = junctions
        
        self.bestPair = None
        self.bestPairType = None
        self.largestRDiff = -1


    #Assign genes and exons to the ASM
    def annotate(self, genome, jxn2gene):
        #Get the set of genes and exons for all junctions 
        genes = dict()
        exons = dict()
        for j in self.junctions:
            for g in jxn2gene[j]:
                genes[g] = 1
                for e in g.exons:
                    exons[e] = 1
        
        self.genes = genes.keys()
        self.exons = list()
        for e in exons:
            for j in self.junctions:
                if e.start==j.end+1 or e.end+1==j.start or (e.start==j.start and e.end==j.end):
                    self.exons.append(e)
                    break
        
        self.genes.sort()
        self.exons.sort()
        self.junctions.sort()
        
            
    #Assign read counts to each junction (and possibly exon) in the ASM.    
    def assignCounts(self, jxnCounts, exnCounts, design):
        self.expDesign = design
        self.jxnCounts = dict()
        self.exnCounts = dict()
        
        #init
        exn2keep = list()
        
        for exn in self.exons:
            if exnCounts and exn in exnCounts:
                self.exnCounts[exn] = exnCounts[exn]
                exn2keep.append(exn)
        self.exons = exn2keep #remove exons that have no count
                
        for jxn in self.junctions:
            self.jxnCounts[jxn] = jxnCounts[jxn]
        
        #Fill in missing items
        for exp in design:
            for cond in design[exp]:
                for exn in self.exons:
                    if exp not in self.exnCounts[exn] or cond not in self.exnCounts[exn][exp]:
                        self.exnCounts[exn].addElement(exp, cond, 0)
                    
                #Assign junction counts
                for jxn in self.junctions:
                    if exp not in self.jxnCounts[jxn] or cond not in self.jxnCounts[jxn][exp]:
                        self.jxnCounts[jxn].addElement(exp, cond, 0)
    
    #Assign read counts but according to a permutation (specific to jSpliceBT)
    def assignShuffledCounts(self, jxnCounts, exnCounts, ref, perm, refCond, permCond):
        self.jxnCounts = dict()
        self.exnCounts = dict()
        
        for exn in self.exons:
            if exnCounts and exn in exnCounts:
                cnt = dict()
                for i in range(len(ref)):
                    e = ref[i]
                    cnt[e] = dict()
                    if e in exnCounts[exn] and refCond in exnCounts[exn][e]:
                        cnt[e][refCond] = exnCounts[exn][e][refCond]
                    else:
                        cnt[e][refCond] = 0
                    p = perm[i]
                    if p in exnCounts[exn] and permCond in exnCounts[exn][p]:
                        cnt[e][permCond] = exnCounts[exn][p][permCond]
                    else:
                        cnt[e][permCond] = 0
                self.exnCounts[exn] = cnt
            else:
                self.exnCounts[exn] = {refCond:0, permCond:0}
                
        for jxn in self.junctions:
            if jxnCounts and jxn in jxnCounts:
                cnt = dict()
                for i in range(len(ref)):
                    e = ref[i]
                    cnt[e] = dict()
                    if e in jxnCounts[jxn] and refCond in jxnCounts[jxn][e]:
                        cnt[e][refCond] = jxnCounts[jxn][e][refCond]
                    else:
                        cnt[e][refCond] = 0
                    p = perm[i]
                    if p in jxnCounts[jxn] and permCond in jxnCounts[jxn][p]:
                        cnt[e][permCond] = jxnCounts[jxn][p][permCond]
                    else:
                        cnt[e][permCond] = 0
                self.jxnCounts[jxn] = cnt
            else:
                self.jxnCounts[jxn] = {refCond:0, permCond:0}
                        

    #Compute the largest relative log2 ratio across all pairs of junctions (or exons) in the ASM.
    def getLargestRDiff(self, rdiffType, countThreshold, inclThreshold, nbExpThreshold, pvalThreshold, rdiffThreshold):
        if rdiffType==ASM.JXN:
            elmnts = self.junctions
            elmntCounts = self.jxnCounts
        elif rdiffType==ASM.EXN:
            elmnts = self.exons
            elmntCounts = self.exnCounts
        else:
            elmnts = list()
        
        
        if len(elmnts)<2:
            return None, None
        
        #Get max reads in all experiments
        maxReads = dict()
        for exp in self.expDesign:
            maxReads[exp] = [0,0]
            for el in elmnts:
                e0 = elmntCounts[el][exp][self.expDesign.cond0]
                e1 = elmntCounts[el][exp][self.expDesign.cond1]
                if rdiffType==ASM.EXN:
                    e0 = e0*1.0/len(el)
                    e1 = e1*1.0/len(el)
        
                if e0>maxReads[exp][0]:
                    maxReads[exp][0] = e0
                if e1>maxReads[exp][1]:
                    maxReads[exp][1] = e1
        
        
        #Select valid elements based on count and inclusion threshold 
        valid = dict()
        for exp in self.expDesign:
            for el in elmnts:
                e0 = elmntCounts[el][exp][self.expDesign.cond0]
                e1 = elmntCounts[el][exp][self.expDesign.cond1]
                ok=False
                if rdiffType==ASM.EXN:
                    ok = (e0>=countThreshold and e0*1.0/len(el)>=inclThreshold*maxReads[exp][0]) or (e1>=countThreshold and e1*1.0/len(el)>=inclThreshold*maxReads[exp][1])
                elif rdiffType==ASM.JXN:
                    ok = (e0>=countThreshold and e0>=inclThreshold*maxReads[exp][0]) or (e1>=countThreshold and e1>=inclThreshold*maxReads[exp][1])
                    
                if ok:
                    if el not in valid:
                        valid[el] = [exp]
                    else:
                        valid[el].append(exp)
                   
        #special case for retained introns
        ri = None
        if rdiffType==ASM.EXN and len(self.junctions)==1:
            for e in elmnts:
                if e==self.junctions[0] and e in valid and len(valid[e])>=nbExpThreshold:
                    ri = e
                    break
            if ri is None: #The retained intron has to be a valid one
                return None,None 
                
        
        #Compute the largest average ratio diff in pairs of elements. (Restricted to overlapping elements if rdiffType=JXN)
        nbElmnts = len(elmnts)
        ratios = dict() #junction -> exps
        pvals = dict() #junction pair -> exps
        rdiffs = dict() #junction pair -> exps
        
        for i in range(nbElmnts-1):
            ei = elmnts[i]
            if ei not in valid or len(valid[ei])<nbExpThreshold:
                continue
            
            for j in range(i+1,nbElmnts):
                ej=elmnts[j]
                if ej not in valid:
                    continue
                
                #special case for retained intron: one element has to be the retained intron
                if ri is not None and ei!=ri and ej!=ri:
                    continue
                
                valExps = list()
                for exp in valid[ei]:
                    if exp in valid[ej]:
                        valExps.append(exp)

                if len(valExps)<nbExpThreshold or (rdiffType==ASM.JXN and not ei.overlaps(ej)):
                    continue

                p = bio.Pair(i,j)
                
                for exp in valExps:
                    e00 = elmntCounts[ei][exp][self.expDesign.cond0]+1
                    e01 = elmntCounts[ei][exp][self.expDesign.cond1]+1
                    e10 = elmntCounts[ej][exp][self.expDesign.cond0]+1
                    e11 = elmntCounts[ej][exp][self.expDesign.cond1]+1
                    
                    r0 = bio.logR(e00, e01)
                    r1 = bio.logR(e10, e11)
                    
                    if math.fabs(r0-r1)>=rdiffThreshold:
                        if r0-r1>0:
                            alt='greater'
                        else:
                            alt='less'
                        _,pval = fisher_exact([[e00,e01],[e10,e11]], alternative=alt) #make a Fisher exact test
                        if pval<=pvalThreshold:
                            if ei not in ratios:
                                ratios[ei] = {exp:r0}
                            else:
                                ratios[ei][exp] = r0
                            if ej not in ratios:
                                ratios[ej] = {exp:r1}
                            else:
                                ratios[ej][exp] = r1
                            
                            if p not in pvals:
                                pvals[p] = dict()
                            pvals[p][exp] = pval
                            
                            if p not in rdiffs:
                                rdiffs[p] = dict()
                            rdiffs[p][exp] = r0-r1
                        
                        
                    
                if p not in rdiffs or len(rdiffs[p])<nbExpThreshold:
                    continue
                
                nbPlus=0
                nbMinus=0
                for r in rdiffs[p].values():
                    if r>0:
                        nbPlus += 1
                    elif r<0:
                        nbMinus += 1
                
                avgdiff = sum(rdiffs[p].values())/len(rdiffs[p])
                absavgdiff = math.fabs(avgdiff)
                if (nbPlus>=nbExpThreshold or nbMinus>=nbExpThreshold) and absavgdiff>self.largestRDiff:
                    self.bestPairType = rdiffType
                    self.largestRDiff = absavgdiff
                    if avgdiff<0:
                        self.bestPair = p.rev()
                    else:
                        self.bestPair = p
        
        #Compute ratio diff for all experiences
        if self.bestPair is None or self.bestPairType!=rdiffType:
            return None, None
        
        if self.bestPair in rdiffs:
            return rdiffs[self.bestPair], pvals[self.bestPair]
        else:
            return rdiffs[self.bestPair.rev()], pvals[self.bestPair.rev()] 
        



    