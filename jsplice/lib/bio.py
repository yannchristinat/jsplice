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

#Biological functions and classes
#For format descriptions, go to https://genome.ucsc.edu/FAQ/FAQformat.html

import math
from expdesign import AbstractDesign

#compute log2(a/b)
def logR(a,b): 
    if a==0:
        if b==0: #both 0
            return float('NaN')
        else: # 0 -> log2(0) = -Inf
            return -float('Inf')
    else:
        if b==0: # log2(a/0) = Inf
            return float('Inf')
        else:
            return math.log(a*1.0/b, 2)
        
def RPKM(cnt, length, totalCnt):
    if length is None:
        return None
    elif cnt==0 or totalCnt==0:
        return 0
    else:
        return cnt*1000.0/(length*totalCnt*1./1000000)
        
#Pair class with comparison and utility functions
class Pair:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def __hash__(self):
        return hash(str(self))
    def __str__(self):
        return str(self.a)+'|'+str(self.b)
    def __eq__(self,p):
        return self.a==p.a and self.b==p.b
    
    def rev(self):
        return Pair(self.b, self.a)


class RegionException(Exception):
    def __init__(self, error):
        super(RegionException, self).__init__('Region: '+error)

#Genomic region class and utility functions
class Region:
    def __init__(self, chrom, start, end, strand):
        if start > end:
            raise RegionException('Start > End')
        
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        
    #Test if the two region overlaps
    def overlaps(self,reg):
        if self.chrom!=reg.chrom or self.strand!=reg.strand: 
            return 0
        if self.start>reg.end or self.end<reg.start:
            return 0
        
        start = self.start
        if self.start<reg.start:
            start = reg.start
        end = self.end
        if self.end>reg.end:
            end = reg.end
        return end-start+1
    
    #Return true if reg is contained in the region
    def contains(self,reg):
        if self.chrom!=reg.chrom or self.strand!=reg.strand: 
            return False
        if self.start<=reg.start and self.end>=reg.end: 
            return True
        return False
    
    def __eq__(self,reg):
        sameChrom=self.chrom==reg.chrom
        if not sameChrom:
            if self.chrom[:3].lower()=='chr' and self.chrom[3:]==reg.chrom:
                sameChrom=True
            elif reg.chrom[:3].lower()=='chr' and reg.chrom[3:]==self.chrom:
                sameChrom=True
        return sameChrom and self.strand==reg.strand and self.start==reg.start and self.end==reg.end
    
    def __hash__(self):
        return hash(str(self))
    
    def __len__(self):
        return self.end-self.start +1
    
    def __str__(self):
        return str(self.chrom)+':'+str(self.start)+'-'+str(self.end)+':'+str(self.strand)
        
    def __lt__(self,reg):
        if self.chrom!=reg.chrom:
            c1 = self.chrom
            c2 = reg.chrom

            # remove chr if present
            if c1[:3]=='chr':
                c1 = c1[3:]
            if c2[:3]=='chr':
                c2 = c2[3:]

            try: #try a numerical comparison
                c1n = int(c1)
                c2n = int(c2)
                return c1n < c2n
            except ValueError: # if the chromosomes are not number, then revert to string comparison
                return c1<c2
        else:
            if self.start<reg.start:
                return True
            elif self.start==reg.start:
                return self.end<reg.end
            else:
                return False

#Gene class. The region of the gene is defined by its exons
#A gene has to be instanciated by one single exon. Subsequent exons are to be added through the 'addExon' method
class Gene(Region):
    def __init__(self, exon, geneId, geneName):
        Region.__init__(self, exon.chrom, exon.start, exon.end, exon.strand)
        self.id = geneId
        self.name = geneName
        self.exons = [exon]
        self.junctions = list()
        self.longestTranscriptLength = None 
        self.rpkm = None
        self.transcripts = list()
        
        
    def addExon(self, exon):
        if exon.start<self.start:
            self.start = exon.start
        if exon.end>self.end:
            self.end=exon.end
        self.exons.append(exon)
        
    def __str__(self):
        return self.id+'|'+self.name+'|'+Region.__str__(self)
    
    def getRegion(self):
        return Region(self.chrom, self.start, self.end, self.strand)
    
    def addRPKM(self, exp, cond, rpkmVal):
        if self.rpkm is None:
            self.rpkm = AbstractDesign()
            
        self.rpkm.addElement(exp,cond,rpkmVal)
        
#Transcript class. The region of the transcript is defined by its exons
#A transcript has to be instanciated by one single exon. Subsequent exons are to be added through the 'addExon' method
class Transcript(Region):
    def __init__(self,exon,trId):
        Region.__init__(self, exon.chrom, exon.start, exon.end, exon.strand)
        self.id = trId
        self.exons = [exon]
    
    def addExon(self, exon):
        if exon.start<self.start:
            self.start = exon.start
        if exon.end>self.end:
            self.end=exon.end
        self.exons.append(exon)
        
            
#Genome class to store a list of genes
class Genome:
    def __init__(self,name):
        self.name = name
        self.genes = list()
        
    #Saves genes and exon position in BED format (sorted). (Designed for jSplice to run coverageBed)
    def saveBED(self, filename):
        # Get all regions
        regs = dict()
        for gene in self.genes:
            regs[gene.getRegion()] = gene.id+'\t0\t'+gene.strand
            for exn in gene.exons:
                regs[exn] = gene.id+'|EXN\t0\t'+exn.strand

        f = open(filename,'w')
        for reg in sorted(regs.keys()):
            f.write(reg.chrom + '\t' + str(reg.start) + '\t' + str(reg.end) + '\t' + regs[reg] + '\n')
        f.close()
        
    def saveGTF(self, filename):
        f = open(filename,'w')
        for gene in self.genes:
            f.write('\t'.join([gene.chrom,'jSplice','CDS',str(gene.start),str(gene.end),'.',gene.strand,'.','gene_id "'+gene.id+'"; transcript_id "'+gene.id+'"; gene_name "'+gene.name+'";'])+'\n')
            exnNb=0
            for exn in sorted(gene.exons):
                exnNb += 1
                f.write('\t'.join([exn.chrom,'jSplice','exon',str(exn.start),str(exn.end),'.',exn.strand,'.','gene_id "'+gene.id+'"; transcript_id "'+gene.id+'"; gene_name "'+gene.name+'"; exon_number '+str(exnNb)+';'])+'\n')
        f.close()
        
        
#Return a genome object with genes and associated exons.
#The file as to in GTF format.  
def processGTFAnnotations(annFile, exonKeyword):
    genes = dict() #gene ID -> gene object
    transcripts = dict() #gene ID -> transcript ID -> list of exons (bio.Region)
    lineNb = 0
    for l in open(annFile,'r'):
        lineNb+=1
        toks = l.strip().split('\t')
        try:
            chrom = toks[0]            
            feature = toks[2]
            start = int(toks[3])
            end = int(toks[4])
            strand = toks[6]
            attrs = toks[8].split(';')
            
            if feature!=exonKeyword or start>end:
                continue
            
            geneId=None
            trId=None
            geneName=''
            for a in attrs:
                a=a.strip()
                if a[:7]=='gene_id': #Mandatory field in GTF files
                    geneId = a[9:-1]
                if a[:13]=='transcript_id': #Mandatory field in GTF files
                    trId = a[15:-1]
                if a[:9]=='gene_name': #Optional field
                    geneName = a[11:-1] 
                
            if geneId is None or trId is None:
                print 'ERROR: Missing gene_id or transcript_id field in GTF file at line '+str(lineNb)
                continue #Mandatory field in GTF files
            
            #Assign/create gene and exon
            exn = Region(chrom,start,end,strand)
            if geneId in genes:
                genes[geneId].addExon(exn)
            else:
                #Create new gene
                gene = Gene(exn,geneId,geneName)
                genes[geneId] = gene
            #Store transcripts
            if geneId not in transcripts:
                transcripts[geneId] = {trId: [exn]}
            elif trId not in transcripts[geneId]:
                transcripts[geneId][trId] = [exn]
            else:
                transcripts[geneId][trId].append(exn)
          
        except:
            continue
        
    #Save genes to genome and compute the longest transcript
    genome = Genome(annFile)
    for geneId in genes:
        gene = genes[geneId]
        #compute longest transcript
        for trId in transcripts[geneId]:
            length=0
            for e in transcripts[geneId][trId]:
                length += len(e)
            if length>gene.longestTranscriptLength or gene.longestTranscriptLength is None:
                gene.longestTranscriptLength = length
        
        genome.genes.append(gene)
    
    if genome.genes:
        return genome
    else:
        print 'ERROR: No genes found in file "'+annFile+'"'
        print 'Check if the file is empty or if it is a valid GTF file!'
        return None

    
