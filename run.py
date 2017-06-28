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

import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import argparse, pickle, warnings
from lib import jsplicelib

warnings.filterwarnings("ignore")

#User parameters
argparser = argparse.ArgumentParser(description='jSplice: a fast method to detect differential alternative splicing events.')
argparser.add_argument('-o','--outdir',help='Output directory', required=True)
argparser.add_argument('-a','--annotation',help='GTF genome annotation file.')
argparser.add_argument('-d','--design', help='Experimental design filename. (Required for the first execution.)')
argparser.add_argument('-t','--type',help='Exon keyword in GTF annotation file (e.g. "exon" in Ensembl GTF files, default)',default='exon')
argparser.add_argument('-e','--exon', help='Flag for exon-based junction files',action='store_true')
argparser.add_argument('-j','--jxnonly',help='Flag for an analysis on junctions only. If set then coverageBed is not performed.', action='store_true')
argparser.add_argument('-b','--nobam',help='Flag for coverageBed. If present, coverageBed is not performed and BED count files are expected instead of BAM files.', action='store_true')
argparser.add_argument('-n','--nbcores',help='Number of CPU cores to use for coverageBed (default: one per BAM file). Used only for coverageBed, hence not relevant if -b or -j is set.',type=int,default=0)
argparser.add_argument('-s','--nostrand',help='Flag for non strand-specific RNA-seq (strand-specific by default). Used only for coverageBed, hence not relevant if -b or -j is set.',action='store_true')
argparser.add_argument('-c','--count',help='Count threshold on junctions and exons (default = 20)',type=int, default=20)
argparser.add_argument('-r','--relfc',help='Relative fold-change threshold (default = 2)',type=float, default=2.0)
argparser.add_argument('-k','--rpkm',help='Gene RPKM threshold (default = 1)',type=float, default=1.0)
argparser.add_argument('-i','--incl',help='Inclusion percentage threshold (default = 0.1)',type=float, default=0.1)
argparser.add_argument('-p','--pvalue',help='Adjusted p-value threshold (default = 0.05)',type=float,default=0.05)
argparser.add_argument('-x','--nbexps',help='Number of experiments that should pass the ratio threshold (all by default)',type=int, default=0)

args = argparser.parse_args()



errmsg=''
if args.relfc<1:
    errmsg += 'ERROR: The relative fold-change threshold (-r,--relfc) should be >= 1\n'
if args.incl<0 or args.incl>1:
    errmsg += 'ERROR: The inclusion percentage threshold (-i,--incl) should be in [0,1]\n'
if args.pvalue<=0 or args.pvalue>1:
    errmsg += 'ERROR: The p-value threshold (-p,--pvalue) should be in ]0,1]\n'
if args.count<1:
    errmsg += 'ERROR: The read count threshold (-c,--count) should be >0\n'
if args.rpkm<0:
    errmsg += 'ERROR: The RPKM threshold (-k,--rpkm) should be >=0\n'
if args.nbcores and args.nbcores<1:
    errmsg += 'ERROR: The number of CPUs (-n) cannot be smaller than 1\n'

#Check outdir and create it if necessary
if os.path.exists(args.outdir):
    if not os.path.isdir(args.outdir):
        errmsg += 'ERROR: '+args.outdir+' is not a directory! Please pick another name.\n'
        sys.exit(0)
else:
    os.makedirs(args.outdir)

#Check if first run
firstrun = not os.path.exists(os.path.join(args.outdir,'jSplice.dat'))
if firstrun and not args.design:
    errmsg += 'ERROR: First run! Please provide an experimental design file. (Not needed for subsequent runs.)\n'

if errmsg:
    print errmsg
    sys.exit(0)

args.outdir = os.path.abspath(args.outdir) #Convert to absolute path

if args.jxnonly:
    print 'jSplice has transformed into Buzz Lightning! (JUNCTION-ONLY mode)'

if not os.path.exists(os.path.join(args.outdir,'jSplice.dat')):
    #run the first step
    genome, expDesign, jxnCounts, exnCounts, totCounts, jxn2gene = jsplicelib.runPreparationStep(args)
else:
    #load jSplice object
    print 'Loading data from previous execution..'
    objs = pickle.load(open(os.path.join(args.outdir,'jSplice.dat'),'rb'))
    genome = objs[0]
    expDesign = objs[1]
    jxnCounts = objs[2]
    exnCounts = objs[3]
    totCounts = objs[4]
    jxn2gene = objs[5]

    if exnCounts is None and not args.jxnonly:
        #Exon counts are missing, hence re-run the first step
        print 'WARNING: The previous run was set to JUNCTION-ONLY. Re-running jSplice from start..'
        genome, expDesign, jxnCounts, exnCounts, totCounts, jxn2gene = jsplicelib.runPreparationStep(args)

if args.nbexps==0:
    args.nbexps = len(expDesign.keys())
elif args.nbexps>len(expDesign.keys()):
    print 'WARNING: The selected number of experiments (-x) is larger than the total number of experiments! (Now set to the maximum.)'
    args.nbexps = len(expDesign.keys())

jsplicelib.runAnalysisStep(genome, expDesign, jxnCounts, exnCounts, totCounts, jxn2gene, args)
