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


# Print the uniquely mapped reads by STAR to a BED6 format
import sys,argparse,os.path

argparser = argparse.ArgumentParser(description='Converts a STAR junction file into a 6-column BED file.')
argparser.add_argument('-m','--multi',help='Flag to sum uniquely mapped junction reads and multi-mapping reads. (Default is to keep only uniquely mapped reads.)',action='store_true')
argparser.add_argument('-f','--starfile',help='STAR junction file',required=True)
argparser.add_argument('-o','--outfile',help='Output file name. (If omitted, the same file name but with a .bed extension will be used.)')

args = argparser.parse_args()

if not os.path.exists(args.starfile):
	print 'ERROR file '+args.starfile+' does not exists!'
	sys.exit()
	
if not args.outfile:
	l=len(args.starfile)
	ext=args.starfile[l-4:]
	if ext=='.csv' or ext=='.tab' or ext=='.txt' or ext=='.tsv':
		args.outfile = args.starfile[:l-4]+'.bed'
	else:
		args.outfile = args.starfile+'.bed'

out=open(args.outfile,'w')
for l in open(args.starfile,'r'):
	toks=l.strip().split('\t')
	cnt=int(toks[6])
	if args.multi:
		cnt += int(toks[7])
	
	strand = '-'
	if toks[3]=='1':
		strand = '+'
		
	out.write("\t".join([toks[0], toks[1], toks[2], 'jxn', str(cnt), strand]) + '\n')
		
out.close()
print 'Output written in '+args.outfile
