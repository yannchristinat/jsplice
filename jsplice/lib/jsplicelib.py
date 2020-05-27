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

import subprocess, Queue, sys, os, pickle, threading, operator, math
import bio, quicksect
from expdesign import AbstractDesign, ExperimentalDesign
from asm import ASM
from bio import Region, logR, RPKM


def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in xrange(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l

#Run coverageBed. Assumes that BAM files are sorted
def runCovBed(q, bedFile, strand_arg):
    #get bedtools version
    cmd = "coverageBed -h 2>&1 >/dev/null | grep Version | awk '{ print $2 }'"
    p=subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    bedtools_version = float('.'.join(out[1:].split('.')[:2]))
    
    while True:
        try:
            exp = q.get(False)
        except Queue.Empty:
            break
        
        if bedtools_version<2.24:
            if strand_arg is None:
                cmd='coverageBed -split -abam "'+exp.bamFile+'" -b "'+bedFile+'" > "'+exp.bedCountFile+'"'
            else:
                cmd='coverageBed '+strand_arg+' -split -abam "'+exp.bamFile+'" -b "'+bedFile+'" > "'+exp.bedCountFile+'"'
        else: #starting on version 2.24, the a and b file are interverted...
            if strand_arg is None:
                cmd='coverageBed -sorted -split -a "'+bedFile+'" -b "'+exp.bamFile+'" > "'+exp.bedCountFile+'"'
            else:
                cmd='coverageBed '+strand_arg+' -sorted -split -a "'+bedFile+'" -b "'+exp.bamFile+'" > "'+exp.bedCountFile+'"'
            
        print cmd
        p=subprocess.Popen(cmd,shell=True) #Not the safest solution but cannot find any other... :(
        p.wait()

# Reads in the junction read counts files
def loadJunctions(expDesign, intronBased):
    jxnCounts = dict() #junction -> AbstractDesign
    for exp in expDesign:
        for cond in expDesign[exp]:
            jxnFile = expDesign[exp][cond].junctionFile
            print 'Loading junction file "'+jxnFile+'"..'

            jxns = dict()
            for l in open(jxnFile,'r'):
                toks = l.strip().split('\t')
                try:
                    chrom=toks[0]
                    start=int(toks[1])
                    end=int(toks[2])
                    cnt=float(toks[4])
                    strand=toks[5]
                    if len(toks)==12: #BED12 -> junctions represented by blocks (TopHat file)
                        blockSizes = toks[10].split(',')
                        start += int(blockSizes[0])
                        end -= int(blockSizes[1])-1

                    #basic check
                    if start>end:
                        print 'ERROR: Invalid junction. Start > end.'
                        continue

                    if not intronBased: #transform to intron-based (historical reasons: the rest of the code is based on intron-based junctions)
                        r = bio.Region(chrom,start+1,end-1,strand)
                    else:
                        r = bio.Region(chrom,start,end,strand)
                    jxns[r]=cnt
                except:
                    continue

            if not jxns:
                print 'ERROR: No elements were found in "'+jxnFile+'".'
                print 'Check if the file is empty or if it is a valid BED file!'

            for j in jxns:
                if j not in jxnCounts:
                    jxnCounts[j] = AbstractDesign()
                jxnCounts[j].addElement(exp, cond, jxns[j])

    return jxnCounts

# Remove junctions without any reads in any sample
def trimZeroJunctions(jxnCounts):
    jxn2del = list()
    for j in jxnCounts:
        valid = False
        for exp in jxnCounts[j]:
            for cond in jxnCounts[j][exp]:
                if jxnCounts[j][exp][cond]>0:
                    valid = True
        if not valid:
            jxn2del.append(j)
    for j in jxn2del:
        del jxnCounts[j]

# Browse the annotation and assign gene(s) to each junction
def assignGenesToJunctions(genome, jxnCounts):
    #create interval tree for all genes
    geneTree = quicksect.IntervalTree()
    for g in genome.genes:
        geneTree.insert(g)

    #Assign genes to junctions
    jxn2gene = dict()
    for j in jxnCounts:
        glist = geneTree.overlaps(j)
        if glist is None:
            jxn2gene[j] = list()
        else:
            jxn2gene[j] = glist
            for g in jxn2gene[j]:
                g.junctions.append(j)

    nogenes=0
    for j in jxn2gene:
        if len(jxn2gene[j])==0:
            nogenes += 1
    print str(nogenes)+' junctions could not be assigned to at least one gene.'
    return jxn2gene

# Create the new genome: Removes exons and genes without any reads and creates new exons through the intersection of junctions and
# annotated exons.
def createTrimmedGenome(genome):
    trimmedGenome = bio.Genome('Trimmed-'+genome.name)

    for gene in genome.genes:
        if len(gene.junctions)==0:
            continue

        trimmedGenome.genes.append(gene) #A gene doesn't need to have exons to be included, just valid junctions

        jxnStarts = dict() #Dummy dictionary
        jxnEnds = dict() #Dummy dictionary
        for j in gene.junctions:
            jxnStarts[j.start] = 1
            jxnEnds[j.end] = 1
        jxnStarts = sorted(jxnStarts)
        jxnEnds = sorted(jxnEnds)

        validExons = dict()
        for exn in gene.exons:
            possibleExnStarts = [exn.start]
            possibleExnEnds = [exn.end]
            for s in jxnStarts:
                if s<=exn.end+1 and s>=exn.start+1 and s-1 not in possibleExnEnds:
                    possibleExnEnds.append(s-1)
                else:
                    break
            for e in jxnEnds:
                if e<=exn.end-1 and e>=exn.start-1 and e+1 not in possibleExnStarts:
                    possibleExnStarts.append(e+1)
                else:
                    break

            #Get all sub-exons
            for start in possibleExnStarts:
                for end in possibleExnEnds:
                    if start>=end: #Invalid exon
                        continue
                    if start==exn.start and end==exn.end: #That is the original exon
                        validExons[exn] = True
                    else:
                        #Create the new exon
                        newExn = bio.Region(exn.chrom, start, end, exn.strand)
                        validExons[newExn] = True

            #Test for retained intron (That is, there exists a junction that is within in the exon)
            for j in gene.junctions:
                if j.start>=exn.start and j.end<=exn.end: #retained intron
                    if j.start==exn.start and j.end==exn.end: #original exon
                        validExons[exn]=True
                    else:
                        #Create new exon
                        newExn = bio.Region(exn.chrom, j.start, j.end, exn.strand)
                        validExons[newExn] = True

        if validExons:
            gene.exons = validExons.keys()


    return trimmedGenome


#Return a dictionary with elements (bio.Region) as keys and read counts as values.
#The file as to be the result of 'coverageBed'.
def getExonScores(countFile):
    genes = dict()
    exons = dict()
    for l in open(countFile,'r'):
        toks = l.strip().split('\t')
        try:
            chrom=toks[0]
            start=int(toks[1])
            end=int(toks[2])
            name = toks[3]
            strand=toks[5]
            cnt=float(toks[6])

            r = Region(chrom,start,end,strand)

            if name.find('|EXN')!=-1:
                exons[r]=cnt
            else:
                genes[r]=cnt
        except:
            continue

    if genes and exons:
        return genes, exons
    else:
        print 'ERROR: No elements were found in "'+countFile+'".'
        print 'Check if the file is empty or if it is a valid BED file!'
        return None, None


# Detect the type of AS event.
# (Note that as we look only at junctions, a true alt-site event cannot be distinguished from an alternative start or end exon.)
def findASMtype(a):
    asmEvent='Unknown ('+str(len(a.junctions))+' junctions)'
    if len(a.junctions)==2:
        j0=a.junctions[0]
        j1=a.junctions[1]
        if j0.start==j1.start:
            if j0.strand=='+':
                asmEvent='Alt 3\' site'
            else:
                asmEvent='Alt 5\' site'
        elif j0.end==j1.end:
            if j0.strand=='+':
                asmEvent='Alt 5\' site'
            else:
                asmEvent='Alt 3\' site'
    elif len(a.junctions)==3:
        srtJxns = sorted(a.junctions)
        j0=srtJxns[0]
        j1=srtJxns[1]
        j2=srtJxns[2]
        if j0.start==j1.start and j1.end==j2.end and j0.end<j2.start:
            asmEvent = 'Cassette exon'
    elif len(a.junctions)==1:
        if a.junctions[0] in a.exons:
            asmEvent = 'Retained intron'
    elif len(a.junctions)==4:
        srtJxns = sorted(a.junctions)
        j0=srtJxns[0]
        j1=srtJxns[1]
        j2=srtJxns[2]
        j3=srtJxns[3]
        if j0.start==j1.start and j2.end==j3.end and j0.end<j2.start and j1.end<j3.start and j2.start<j1.end:
            asmEvent = 'Mut. excl. exons'
    if asmEvent[:7]=='Unknown' and len(a.junctions)>2:
        sameLeft=True
        sameRight=True
        for j in a.junctions:
            if j.start!=a.junctions[0].start:
                sameLeft = False
            if j.end!=a.junctions[0].end:
                sameRight = False

        if sameLeft:
            if a.junctions[0].strand=='+':
                asmEvent = 'Mult. alt. 3\' site'
            else:
                asmEvent = 'Mult. alt. 5\' site'
        if sameRight:
            if a.junctions[0].strand=='+':
                asmEvent = 'Mult. alt. 5\' site'
            else:
                asmEvent = 'Mult. alt. 3\' site'

    return asmEvent

# Write results into a single HTML file. (Requires an internet connection to retrieve the jQuery script.)
def printHtml(htmlfilename, srtAvgRDiffs, rdiffs, pvals, exonCnts, junctionCnts, totCounts, args, expdesign):
    f = open(htmlfilename,'w')

    exps = expdesign.keys()
    c0 = expdesign.cond0
    c1 = expdesign.cond1

    expStr = ''
    for e in exps:
        expStr += ', '+str(e)
    expStr = expStr[2:]
    nbExps = len(exps)


    #
    # Write header, scripts and stuff
    #
    f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    f.write('<html xmlns="http://www.w3.org/1999/xhtml" >\n')
    f.write('<head>\n')
    f.write('<title>jSplice results</title>\n')
    f.write('<style type="text/css">\n')
    f.write('        body { font-family:Arial, Helvetica, Sans-Serif; font-size:0.8em;}\n')
    f.write('        #report { border-collapse:collapse;}\n')
    f.write('        #report h4 { margin:0px; padding:0px;}\n')
    f.write('        #report img { float:right;}\n')
    f.write('        #report ul { margin:10px 0 10px 40px; padding:0px;}\n')
    f.write('        #report th { background:#624d41 none repeat-x scroll center left; color:#fff; padding:5px 15px; text-align:left;}\n')
    f.write('        #report td { background:#f6f3e8 none repeat-x scroll center left; padding:2px 15px; }\n')
    f.write('        #report tr.click td { background:#ded2ba none repeat-x scroll center left; cursor:pointer; }\n')
    #f.write('        #report div.arrow { background:transparent url(arrows.png) no-repeat scroll 0px -16px; width:16px; height:16px; display:block;}\n')
    #f.write('        #report div.up { background-position:0px 0px;}\n')
    f.write('    </style>\n')
    f.write('    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js" type="text/javascript"></script>\n')
    f.write('    <script type="text/javascript">\n')
    f.write('        $(document).ready(function(){\n')
    f.write('            $("#report tr.details").hide();\n')
    f.write('            $("#report tr.click").click(function(){\n')
    f.write('                $(this).next("tr").toggle();\n')
    #f.write('                $(this).find(".arrow").toggleClass("up");\n')
    f.write('            });\n')
    f.write('            //$("#report").jExpand();\n')
    f.write('        });\n')
    f.write('    </script>\n')
    f.write('</head>\n')
    f.write('<body>\n')
    f.write('<h1>jSplice</h1>\n')


    #
    # Write parameters
    #
    f.write('<h3>Parameters</h3>\n')
    f.write('<table>\n')
    if args.jxnonly:
        f.write('<tr><td>Results based on:</td><td>Junctions only</td></tr>\n')
    else:
        f.write('<tr><td>Results based on:</td><td>Junctions + Exons</td></tr>\n')
    f.write('<tr><td>Count threshold:</td><td>'+str(args.count)+'</td></tr>\n')
    f.write('<tr><td>Relative fold-change threshold:</td><td>'+str(args.relfc)+'</td></tr>\n')
    f.write('<tr><td>Gene RPKM threshold:</td><td>'+str(args.rpkm)+'</td></tr>\n')
    f.write('<tr><td>Inclusion percentage threshold:</td><td>'+str(args.incl*100)+'%</td></tr>\n')
    f.write('<tr><td>P-value threshold:</td><td>'+str(args.pvalue)+'</td></tr>\n')
    f.write('<tr><td>Number of experiments threshold:</td><td>'+str(args.nbexps)+'/'+str(len(expdesign.keys()))+'</td></tr>\n')
    f.write('<tr><td>Results directory:</td><td>'+str(args.outdir)+'</td></tr>\n')
    f.write('</table>\n')
    f.write('<h3>Results</h3>\n')
    f.write('<p style="font-style:italic">Note: exons or junctions highlighted in red were used to compute the ratio difference.</p>\n')


    #
    # Write ASMs
    #
    f.write('<table id="report">\n')
    f.write('<tr>\n')
    f.write('    <th>Gene name</th>\n')
    f.write('    <th>ASM type</th>\n')
    f.write('    <th>Average relative log2('+c0+'/'+c1+')</th>\n')
    f.write('    <th colspan='+str(nbExps)+'><div align=center>Largest relative log2('+c0+'/'+c1+')<br>('+expStr+')</div></th>\n')
    f.write('    <th colspan='+str(nbExps)+'><div align=center>Adjusted p-values (FDR)<br>('+expStr+')</div></th>\n')
    #f.write('    <th></th>\n')
    f.write('</tr>\n')

    for p in srtAvgRDiffs:
        a=p[0]
        avgRDiff=p[1]

        #detect the type of AS event
        asmEvent=findASMtype(a)

        #write main info
        f.write('<tr class="click">\n')

        gnstr = ''
        for g in a.genes:
            if g.name!="":
                gnstr += '<br>'+g.name+'|'+g.id
            else:
                gnstr += '<br>'+g.id
        if gnstr=='':
            gnstr = 'Unknown'
        else:
            gnstr = gnstr[4:]

        f.write('<td>'+gnstr+'</td><td>'+asmEvent+'</td><td>'+str("{0:.3f}".format(avgRDiff))+'</td>\n')
        for e in exps:
            if e in rdiffs[a]:
                f.write('<td>'+str("{0:.3f}".format(rdiffs[a][e]))+'</td>')
            else:
                f.write('<td>N/A</td>')
        f.write('\n')
        for e in exps:
            if a in pvals and e in pvals[a]:
                f.write('<td>'+str("{0:.5f}".format(pvals[a][e]))+'</td>')
            else:
                f.write('<td>N/A</td>')
        f.write('\n')
        #f.write('<td><div class="arrow">&nbsp;</div></td></tr>\n')

        #write detailed info
        f.write('<tr class="details"><td colspan='+str(nbExps*3+2)+'>\n')

        detailTableHeader = '<th colspan='+str(nbExps)+'><div align=center>log2('+c0+'/'+c1+')<br>('+expStr+')</div></th>\n'
        detailTableHeader += '<th colspan='+str(nbExps)+'><div align=center>'+c0+' (RPKM)<br>('+expStr+')</div></th>\n'
        detailTableHeader += '<th colspan='+str(nbExps)+'><div align=center>'+c1+' (RPKM)<br>('+expStr+')</div></th></tr>\n'

        #write gene expression
        f.write('<table>\n')
        f.write('<tr><th style="background:transparent"><h3 style="color:#000">Gene expression</h3></th>'+detailTableHeader+'\n')

        for g in a.genes:
            gidStr=''
            if g.name!="":
                gidStr += g.name+'|'
            if g.id[:4]=='ENSG': #Ensembl gene ID
                gidStr += '<a href="http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g='+g.id+'" target="_blank">'+g.id+'</a>'
            elif g.id[:4]=='ENSMUSG': #Ensembl gene ID
                gidStr += '<a href="http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g='+g.id+'" target="_blank">'+g.id+'</a>'
            else:
                gidStr += g.id

            f.write('<tr><td>'+gidStr+'</td>\n')
            for e in exps:
                if g.rpkm is not None and e in g.rpkm and c0 in g.rpkm[e] and c1 in g.rpkm[e]:
                    f.write('<td>'+str("{0:.3f}".format(logR(g.rpkm[e][c0],g.rpkm[e][c1])))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if g.rpkm is not None and e in g.rpkm and c0 in g.rpkm[e] :
                    f.write('<td>'+str("{0:.1f}".format(g.rpkm[e][c0]))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if g.rpkm is not None and e in g.rpkm and c1 in g.rpkm[e] :
                    f.write('<td>'+str("{0:.1f}".format(g.rpkm[e][c1]))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('</tr>\n')


        #write exon expression
        f.write('<tr><td colspan='+str(nbExps*3+1)+'>&nbsp;</td></tr>\n')
        f.write('<tr><th style="background:transparent"><h3 style="color:#000">Exon expression</h3></th>'+detailTableHeader+'\n')
        for exn in a.exons:
            style=''
            if a.bestPairType==ASM.EXN and (a.exons[a.bestPair.a]==exn or a.exons[a.bestPair.b]==exn):
                style = ' style="color:#e50005"'

            if exn not in exonCnts:
                continue
            f.write('<tr'+style+'><td>'+str(exn)+'</td>\n')
            for e in exps:
                if e in exonCnts[exn] and c0 in exonCnts[exn][e] and c1 in exonCnts[exn][e]:
                    f.write('<td>'+str("{0:.3f}".format(logR(exonCnts[exn][e][c0],exonCnts[exn][e][c1])))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if e in exonCnts[exn] and c0 in exonCnts[exn][e] :
                    rpkm = RPKM(exonCnts[exn][e][c0], len(exn), totCounts[e][c0])
                    f.write('<td>'+str("{0:.1f}".format(rpkm))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if e in exonCnts[exn] and c1 in exonCnts[exn][e] :
                    rpkm = RPKM(exonCnts[exn][e][c1], len(exn), totCounts[e][c1])
                    f.write('<td>'+str("{0:.1f}".format(rpkm))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('</tr>\n')
        f.write('<tr><td colspan='+str(nbExps*3+1)+'>&nbsp;</td></tr>\n')

        #write junction expression
        f.write('<tr><th style="background:transparent"><h3 style="color:#000">Junction expression</h3></th>'+detailTableHeader.replace('RPKM','raw reads')+'\n')
        for jxn in a.junctions:
            style=''
            if a.bestPairType==ASM.JXN and (a.junctions[a.bestPair.a]==jxn or a.junctions[a.bestPair.b]==jxn):
                style = ' style="color:#e50005"'

            if jxn not in junctionCnts:
                continue
            f.write('<tr'+style+'><td>'+str(jxn)+'</td>\n')
            for e in exps:
                if e in junctionCnts[jxn] and c0 in junctionCnts[jxn][e] and c1 in junctionCnts[jxn][e]:
                    f.write('<td>'+str("{0:.3f}".format(logR(junctionCnts[jxn][e][c0],junctionCnts[jxn][e][c1])))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if e in junctionCnts[jxn] and c0 in junctionCnts[jxn][e] :
                    f.write('<td>'+str(int(junctionCnts[jxn][e][c0]))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('\n')
            for e in exps:
                if e in junctionCnts[jxn] and c1 in junctionCnts[jxn][e] :
                    f.write('<td>'+str(int(junctionCnts[jxn][e][c1]))+'</td>')
                else:
                    f.write('<td>N/A</td>')
            f.write('</tr>\n')

        f.write('<tr><td colspan='+str(nbExps*3+1)+'>&nbsp;</td></tr>\n')
        f.write('</table>\n')
        f.write('</td></tr>\n')

    # Write closing tags
    f.write('</table><p>&nbsp;</p><p>If you use this software in your research, please cite: <br/>Christinat Y., Pawlowski R., and Krek W. jSplice: a high-performance method for accurate prediction of alternative splicing events and its application to large-scale renal cancer transcriptome data. Bioinformatics, 2016, 1-9</p><p>Copyright Yann Christinat 2015</p></body>\n')

# Write results into a text file (read counts are not included).
def printText(filename, srtAvgRDiffs, rdiffs, pval):
    outfile = open(filename,'w')
    outfile.write('Gene_name|ID\tASM_type\tAvg_relFC\tLargest_relFCs\tadjPvalues\tASM_junctions\tASM_exons\n')
    for p in srtAvgRDiffs:
        a=p[0]
        avgRDiff = p[1]


        pvalStr=''
        for exp in pval[a]:
            pvalStr += ','+exp+':'+str(pval[a][exp])

        rdiffStr=''
        for exp in rdiffs[a]:
            rdiffStr += ','+exp+':'+str(rdiffs[a][exp])

        jxnStr=''
        for j in a.junctions:
            if a.bestPairType==ASM.EXN and (a.exons[a.bestPair.a]==j or a.exons[a.bestPair.b]==j):
                jxnStr+=',*'+str(j)
            else:
                jxnStr+=','+str(j)

        exnStr=''
        for e in a.exons:
            if a.bestPairType==ASM.EXN and (a.exons[a.bestPair.a]==e or a.exons[a.bestPair.b]==e):
                exnStr+=',*'+str(e) #Add a '*' if the exon was used to compute the
            else:
                exnStr+=','+str(e)

        asmEvent = findASMtype(a)

        gnstr = ''
        for g in a.genes:
            if g.name!="":
                gnstr += ','+g.name+'|'+g.id
            else:
                gnstr += ','+g.id

        outfile.write(gnstr[1:]+'\t'+asmEvent+'\t'+str(avgRDiff)+'\t'+rdiffStr[1:]+'\t'+pvalStr[1:]+'\t'+jxnStr[1:]+'\t'+exnStr[1:]+'\n')
    outfile.close()

# Compute RPKM values for each gene. (Based on the longest transcript of the gene)
def setRPKM(genome, geneCounts, totCount):
    for exp in totCount:
        for cond in totCount[exp]:
            if totCount[exp][cond]==0:
                continue

            for gene in genome.genes:
                genereg = gene.getRegion()
                if genereg not in geneCounts or exp not in geneCounts[genereg] or cond not in geneCounts[genereg][exp]:
                    gene.addRPKM(exp,cond,0)
                else:
                    cnt = geneCounts[genereg][exp][cond]
                    rpkm = RPKM(cnt, gene.longestTranscriptLength, totCount[exp][cond])
                    gene.addRPKM(exp, cond, rpkm)


#Trim junctions
#To be valid a junction has to be above the threshold in one condition of each experiment (or N if 'nbexps' is set)
def trimGenome(genome, jxn2gene, jxnCounts, exnCounts, countThr, rpkmThr, nbexpThr):
    #Trim junctions
    jxn2del = list()
    for j in jxnCounts:
        nbValExps = 0
        if j not in jxnCounts:
            continue
        for exp in jxnCounts[j]:
            valid = False
            for cond in jxnCounts[j][exp]:
                if jxnCounts[j][exp][cond]>=countThr:
                    valid = True
            if valid:
                nbValExps += 1
        if nbValExps<nbexpThr:
            jxn2del.append(j)

    for j in jxn2del:
        for g in jxn2gene[j]:
            g.junctions.remove(j)
        del jxn2gene[j]

    #Trim genes and exons
    gene2keep = quicksect.IntervalTree()
    newgenes = list()
    for gene in genome.genes:
        #First test if there remains any valid junction
        if len(gene.junctions)==0:
            continue

        #Test RPKM (only if possible, if in junction-only mode, then no RPKM or exon count can be used)
        if gene.rpkm is None or exnCounts is None:
            gene.exons = list()
            newgenes.append(gene)
            gene2keep.insert(gene)
        else:
            nbValExps = 0
            for exp in gene.rpkm:
                valid = True
                for cond in gene.rpkm[exp]:
                    if gene.rpkm[exp][cond]<rpkmThr:
                        valid = False
                if valid:
                    nbValExps += 1
            if nbValExps<nbexpThr:
                continue

            #keep gene
            newgenes.append(gene)
            gene2keep.insert(gene)

            #Trim exons
            exn2keep = list()
            for e in gene.exons:
                #Test if there still exists a junction to validate the exon
                valid = False
                for j in gene.junctions:
                    if j.start==e.end+1 or j.end+1==e.start or j==e:
                        valid = True
                        break
                if not valid:
                    continue

                #Test if the read count is higher than the threshold
                if e in exnCounts:
                    nbValExps = 0
                    for exp in exnCounts[e]:
                        for cond in exnCounts[e][exp]:
                            if exnCounts[e][exp][cond]>=countThr:
                                nbValExps += 1
                                break #there exists at least one condition where the threshold is met.
                    valid = nbValExps>=nbexpThr

                if valid:
                    exn2keep.append(e)

            #Update gene
            gene.exons = exn2keep

    #Remove genes from genome
    genome.genes = newgenes


    #Clean jxn2gene
    for j in jxn2gene:
        valgenes=list()
        for g in jxn2gene[j]:
            if gene2keep.contains(g):
                valgenes.append(g)
        jxn2gene[j] = valgenes


def getChromosomeScheme(chromList):
    scheme = None
    if 'chr1' in chromList:
        scheme = 'UCSC'
    elif '1' in chromList:
        scheme = 'Ensembl'
    return scheme

def correctAnnotations(junctions, genome):
    #Get chromosome lists
    jxnChromList = list()
    for j in junctions:
        if j.chrom not in jxnChromList:
            jxnChromList.append(j.chrom)

    annChromList = list()
    for g in genome.genes:
        if g.chrom not in annChromList:
            annChromList.append(g.chrom)

    #Get overlap in chromosomes
    unmatchedChrom = list()
    for chrom in jxnChromList:
        if chrom not in annChromList:
            unmatchedChrom.append(chrom)
    
    print 'Unmatched chrom. in annotation: '+str(unmatchedChrom)
    
    if len(unmatchedChrom)>0:
        print 'WARNING: Chromosome names in annotation file are different from the junctions\'.'

        jxnScheme = getChromosomeScheme(jxnChromList)
        print 'Junctions chrom. scheme: '+jxnScheme
        annScheme = getChromosomeScheme(annChromList)
        print 'Annotation chrom. scheme: '+annScheme
        
        if jxnScheme is not None and annScheme is not None and jxnScheme!=annScheme:
            print 'jSplice will attempt to correct the annotations by adding or removing "chr" in front of the chromosome name.'

            if jxnScheme=='UCSC' and annScheme=='Ensembl': #Add chr to all chromosomes in annotations
                for g in genome.genes:
                    g.chrom = 'chr'+g.chrom
                    for e in g.exons:
                        e.chrom = 'chr'+e.chrom
            elif jxnScheme=='Ensembl' and annScheme=='UCSC': #Remove chr to all chromosomes in annotations
                for g in genome.genes:
                    if len(g.chrom)>3 and g.chrom[:3].lower()=='chr':
                        g.chrom = g.chrom[3:]
                    for e in g.exons:
                        if len(e.chrom)>3 and e.chrom[:3].lower()=='chr':
                            e.chrom = e.chrom[3:]
        else:
            print 'Unknown naming scheme (nor Ensembl, nor UCSC). jSplice will run without annotations (and thus in junction-only mode).'
            return False

    return True

# Reads in the junction and exon read counts files
def loadCounts(expDesign):
    exonCounts = dict() # exon -> AbstractDesign
    geneCounts = dict() # gene -> AbstractDesign
    totalCounts = AbstractDesign() #sum of reads across all genes. Note: reads not mapped to genes are discarded and reads contained in several genes are overcounted.
    for exp in expDesign:
        for cond in expDesign[exp]:
            countFile = expDesign[exp][cond].bedCountFile
            if countFile:
                totalCounts.addElement(exp, cond, 0)

                print 'Loading gene/exon CoverageBED file "'+countFile+'"..'

                geneCnts = dict()
                exnCnts = dict()
                for l in open(countFile,'r'):
                    toks = l.strip().split('\t')
                    try:
                        chrom=toks[0]
                        start=int(toks[1])
                        end=int(toks[2])
                        name = toks[3]
                        strand=toks[5]
                        cnt=float(toks[6])

                        r = Region(chrom,start,end,strand)

                        if name.find('|EXN')!=-1:
                            exnCnts[r]=cnt
                        else:
                            geneCnts[r]=cnt
                    except:
                        continue

                if not geneCnts and not exnCnts:
                    print 'ERROR: No elements were found in "'+countFile+'".'
                    print 'Check if the file is empty or if it is a valid BED file!'


                for e in exnCnts:
                    if e not in exonCounts:
                        exonCounts[e] = AbstractDesign()
                    exonCounts[e].addElement(exp,cond,exnCnts[e])

                for g in geneCnts:
                    if g not in geneCounts:
                        geneCounts[g] = AbstractDesign()
                    geneCounts[g].addElement(exp,cond,geneCnts[g])

                    totalCounts[exp][cond] += geneCnts[g]

    #Sanity check on total counts
    for exp in totalCounts:
        for cond in totalCounts[exp]:
            if totalCounts[exp][cond] == 0:
                print 'WARNING! Experiment '+exp+': '+cond+' has zero reads!'

    return geneCounts,exonCounts,totalCounts

# Find all ASMs in the genome
def findASMs(genome, jxn2gene):
    asms = list() #list of ASMs

    #create a junction tree
    jTree = dict()
    for j in jxn2gene:
        if j.chrom not in jTree:
            jTree[j.chrom] = {j.strand:[j]}
        elif j.strand not in jTree[j.chrom]:
            jTree[j.chrom][j.strand] = [j]
        else:
            jTree[j.chrom][j.strand].append(j)
    #jTree = GTree()
    #for j in jxn2gene:
    #    jTree.add(j)

    #find overlapping junctions
    for chrom in jTree:
        for strand in jTree[chrom]:
            junctions = sorted(jTree[chrom][strand])
            nbJxns = len(junctions)

            asmIds = [0]*nbJxns
            nextASMid = 1
            for a in range(nbJxns-1):
                jA = junctions[a]
                for b in range(a+1, nbJxns):
                    jB = junctions[b]
                    if jB.start>jA.end:
                        break
                    elif jA.overlaps(jB):
                        if asmIds[a]==0:
                            if asmIds[b]==0:
                                asmIds[a]=nextASMid
                                asmIds[b]=nextASMid
                                nextASMid += 1
                            else:
                                asmIds[a]=asmIds[b]
                        else:
                            if asmIds[b]==0:
                                asmIds[b]=asmIds[a]
                            else:
                                oldId = asmIds[a]
                                newId = asmIds[b]
                                for k in range(nbJxns):
                                    if asmIds[k]==oldId:
                                        asmIds[k]=newId

            #Create ASMs
            asmDict = dict()
            for i in range(nbJxns):
                aid = asmIds[i]
                if aid==0:
                    continue

                if aid not in asmDict:
                    asmDict[aid] = list()
                asmDict[aid].append(junctions[i])

            for aid in asmDict:
                a = ASM(asmDict[aid])
                a.annotate(genome, jxn2gene)
                asms.append(a)

            # Find retained introns
            for i in range(nbJxns):
                if asmIds[i]==0: #Only then test for RI
                    j = junctions[i]
                    exons = dict() #contains all exons from all genes associated to this junction
                    for g in jxn2gene[j]:
                        for e in g.exons:
                            exons[e] = 1

                    up=list()
                    ri=None
                    down=list()
                    for e in exons:
                        if e.start==j.end+1:
                            down.append(e)
                        elif e.end+1==j.start:
                            up.append(e)
                        elif e.start==j.start and e.end==j.end:
                            ri = e

                    if ri and up and down:
                        a = ASM([j])

                        #fill in exons and genes
                        a.exons = up + [ri] + down
                        a.genes = jxn2gene[j]

                        asms.append(a)

    return asms



def runPreparationStep(args):
    #Create ExpDesign object and add count file names
    flag = 0
    if args.nobam:
        flag = ExperimentalDesign.NOBAM
    if args.jxnonly:
        flag = ExperimentalDesign.JXNONLY
        
    if not args.design:
        print 'ERROR: No experimental design provided!'
        sys.exit(0)

    expDesign = ExperimentalDesign(args.design, flag)
    if not expDesign.valid:
        print 'ERROR: Experimental design in not valid!'
        sys.exit(0)

    if not args.nobam and not args.jxnonly: #if coverageBed will be run, then create bed count files
        for exp in expDesign:
            for cond in expDesign[exp]:
                expDesign[exp][cond].bedCountFile = os.path.join(args.outdir,exp+'-'+cond+'.counts.bed')
    if args.jxnonly:
        for exp in expDesign:
            for cond in expDesign[exp]:
                expDesign[exp][cond].bedCountFile = None
                expDesign[exp][cond].bamFile = None

    #Load junction files
    jxnCounts = loadJunctions(expDesign, not args.exon)
    trimZeroJunctions(jxnCounts) #In case some junctions have zero counts in all experiments

    if len(jxnCounts)==0:
        print 'ERROR: No junctions in files or all had zero reads!'
        sys.exit(0)
    print str(len(jxnCounts))+' valid junctions across experiments'



    trimmedGenome = None
    jxn2gene = None
    if not args.annotation:
        jxn2gene = dict()
        for j in jxnCounts:
            jxn2gene[j] = list() #fill in jxn2gene
    else:
        #Load GTF annotation
        print 'Loading GTF genome annotation file "'+args.annotation+'"..'
        genomeAnn = bio.processGTFAnnotations(args.annotation, args.type)
        print str(len(genomeAnn.genes))+' genes in GTF file.'

        #Correct annotation
        corr = correctAnnotations(jxnCounts.keys(), genomeAnn)
        if not corr:
            args.jxnonly = True

        else:
            #Assign junctions to genes
            print 'Assigning junctions to genes..'
            jxn2gene = assignGenesToJunctions(genomeAnn, jxnCounts)

            #Create trimmed genome
            print 'Creating new genome..'
            trimmedGenome = createTrimmedGenome(genomeAnn)
            print str(len(trimmedGenome.genes))+' gene in new genome.'
            if len(trimmedGenome.genes)==0:
                print 'ERROR: No genes left in new genome! Check your junction and GTF files.'
                sys.exit()

            genomeFile = os.path.join(args.outdir,'jSplice_genome.bed')
            trimmedGenome.saveBED(genomeFile)
            #trimmedGenome.saveGTF(os.path.join(args.outdir,'jSplice_genome.gtf')) #Experimental!!!!

            print 'New genome is saved in "'+os.path.join(args.outdir,'jSplice_genome.bed')+'"'

            #Running coverageBed (if needed)
            if not args.jxnonly and not args.nobam:
                strand = None
                if args.samestrand:
                    strand = '-s'
                elif args.diffstrand:
                    strand = '-S'
                
                if args.nbcores==0:
                    args.nbcores = len(expDesign)
                print 'Running coverageBed on '+str(args.nbcores)+' threads.. (might take a long time)'

                q = Queue.Queue()
                for exp in expDesign:
                    for cond in expDesign[exp]:
                        e = expDesign[exp][cond]
                        q.put(e)

                tList = list()
                for _ in range(args.nbcores):
                    t = threading.Thread(target=runCovBed, args=(q,genomeFile,strand))
                    t.daemon=True #Thread dies if the main process dies
                    tList.append(t)
                    t.start()

                for t in tList:
                    t.join()

                print 'Done'

    #Load the exon counts if not in JUNCTION ONLY mode
    geneCounts = None
    exnCounts = None
    totCount = None
    if not args.jxnonly and args.annotation:
        geneCounts,exnCounts,totCount = loadCounts(expDesign)
        setRPKM(trimmedGenome, geneCounts, totCount) #Compute RPKM for each gene


    #Saving jSplice object: expDesign, junction counts, new genome
    print 'Saving jSplice objects in "'+os.path.join(args.outdir,'jSplice.dat')+'"'
    pickle.dump([trimmedGenome, expDesign, jxnCounts, exnCounts, totCount, jxn2gene], open(os.path.join(args.outdir,'jSplice.dat'),'wb'), pickle.HIGHEST_PROTOCOL)

    return trimmedGenome, expDesign, jxnCounts, exnCounts, totCount, jxn2gene


def runAnalysisStep(genome, expDesign, junctionCounts, exonCounts, totCount, jxn2gene, args):


    #Trim genome to new count threshold and RPKM threshold
    if genome is not None:
        print 'Trimming genome..'
        trimGenome(genome, jxn2gene, junctionCounts, exonCounts, args.count, args.rpkm, args.nbexps)
        print str(len(genome.genes))+' genes left in trimmed genome'



    # Detect ASMs
    print 'Detecting ASMs..'
    asms = findASMs(genome, jxn2gene)
    print str(len(asms))+' ASMs detected.'

    #Statistics:
    ri=0
    geneless=0
    singlegene=0
    multigene=0
    for a in asms:
        if len(a.junctions)==1:
            ri +=1
        if len(a.genes)==0:
            geneless += 1
        elif len(a.genes)==1:
            singlegene += 1
        elif len(a.genes)>1:
            multigene += 1
    print 'Geneless: '+str(geneless)+', Single-gene: '+str(singlegene)+' Multigene: '+str(multigene)
    print 'Retained introns: '+str(ri)



    # Analyse ASMs
    print 'Analyzing ASMs..'
    rdiffs = dict()
    avgRDiffs = dict()
    allpvals = list()
    pvals=dict()
    for a in asms:
        a.assignCounts(junctionCounts, exonCounts, expDesign)
        rdiffExn = None
        if not args.jxnonly:
            rdiffExn, pvalExn = a.getLargestRDiff(ASM.EXN, args.count, args.incl, args.nbexps ,args.pvalue, math.log(args.relfc,2))
        rdiffJxn, pvalJxn = a.getLargestRDiff(ASM.JXN, args.count, args.incl, args.nbexps, args.pvalue, math.log(args.relfc,2))

        if rdiffJxn:
            rdiffs[a] = rdiffJxn
            avgRDiffs[a] = a.largestRDiff
            pvals[a]=pvalJxn
            allpvals += pvalJxn.values()
        elif not args.jxnonly and rdiffExn:
            rdiffs[a] = rdiffExn
            avgRDiffs[a] = a.largestRDiff
            pvals[a]=pvalExn
            allpvals += pvalExn.values()

    #Compute FDR and trim ASMs
    fdrVals = FDR(allpvals)
    fdrs=dict()
    for i in range(len(allpvals)):
        fdrs[allpvals[i]] = fdrVals[i]
    asmFDRs=dict()
    for a in pvals:
        fdrValList=dict()
        for exp in pvals[a]:
            v = pvals[a][exp]
            fdrValList[exp] = fdrs[v]
        asmFDRs[a] = fdrValList

    #Print summary file
    srtAvgRDiffs = sorted(avgRDiffs.iteritems(), key=operator.itemgetter(1), reverse=True)

    printHtml(os.path.join(args.outdir,'jSplice_results.html'), srtAvgRDiffs, rdiffs, asmFDRs, exonCounts, junctionCounts, totCount, args, expDesign)
    printText(os.path.join(args.outdir,'jSplice_results.txt'),srtAvgRDiffs, rdiffs, asmFDRs)

    print str(len(avgRDiffs))+' ASMs selected. (See file "'+os.path.join(args.outdir,'jSplice_results.html')+'" for details.)'








