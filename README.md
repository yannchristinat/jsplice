jSplice 1.0.3
-------------

jSplice is a fast method to detect differential alternative splicing events from RNA-seq experiments. (see publication for details).

jSplice is freely available under the GPL license at http://www.mhs.biol.ethz.ch/research/krek/jsplice

Please use the following citation if you use jSplice for your research:

Yann Christinat, Rafał Pawłowski, Wilhelm Krek; jSplice: a high-performance method for accurate prediction of alternative splicing events and its application to large-scale renal cancer transcriptome data. Bioinformatics 2016; 32 (14): 2111-2119.

What's new in version 1.0.3
---------------------------
- Fixed the issue with version 2.24.0+ of bedtools (where they changed the command line format)
- Proper packaging of jSplice and system-wide commandline tool
- Extended the strandedness option for coverageBed to include -s and -S as defined by coverageBed

Requirements
------------
- Python 2.7 with numpy/scipy modules
- coverageBed from BEDtools package


Installation
------------
sudo python setup.py install

If you don't have administrator rights, then append "--user" to the command. Also make sure that ~/.local/bin is in your path (or move the files somewhere else).


Running jSplice
---------------

jsplice [-h] [-t KEYWORD] [-e] [-b] [-n N] [-s] [-S] [-j] [-x NBEXPS] [-c COUNT] [-r RDIFF] [-k RPKM] [-i INCL] [-p PVAL] [–a ANNOTATION] –d DESIGN –o OUTDIR

Upon the first run, jSplice saves the results into an object file (“OUTDIR/jSplice.dat”). If one wants to re-run jSplice with different parameters or thresholds, the object is loaded to lower computational costs.

Note: if the jSplice folder is in your PYTHONPATH then you can call it as a module: python –m jsplice.run


Required parameters
-------------------

-o, --outdir: Output directory. Any file already present in this directory will be overwritten. Ideally, the directory should identify your experiment. After execution, the directory will contain a jSplice object, a BED file of the trimmed genome, and all files generated by coverageBed (if any). Note that the directory will be used by runJSplice.py later on.

-d, --design: Path to the experimental design file. The text file should contain one line per experiment and four columns: experiment, condition, junction file, and BAM file. If –b is set, then a BED file is expected in the fourth column instead of the BAM file. If –j is set, then the fourth column can be skipped. 
Important: file paths in the experimental design file can be either absolute or relative but cannot contain “~”.
This parameter is only required for the first run.


Optional parameters - controlling the thresholds
------------------------------------------------

-c, --count: Read count threshold. Any junction or exon that does not reach the threshold in at least one condition for all experiment is discarded. Default is 20.

-r, --relfc: Relative fold-change threshold. The ratio difference of two element of an ASM is computed as the difference between the fold-change of condition A to B. The pair of elements that maximizes the fold-change difference (or relative fold-change) across all experiments is selected as representative of the ASM. The average log2 ratio difference of the representative pair is used to sort ASMs. If its value does not meet the threshold then the ASM is discarded. Note that a fold-change ratio is expected, not its log2 transform. 
(2 by default)

-k, --rpkm: RPKM threshold for gene expression. This parameter allows a simple filtering of non-expressed genes. To be retained, a gene has to be above the threshold in all experiments and all conditions. (1 by default)

-i, --incl: Inclusion percentage threshold (between 0 and 1). To be selected as a part of the representative pair (and thus be used in the computation of the largest average relative fold-change), an element has to be strongly expressed in at least one condition. (0.1 by default, which corresponds to a 10% threshold)


-p, --pvalue: Exact Fisher’s test p-value threshold. (0.05 by default)

-x, --nbexps: Number of experiments to consider for all thresholds. For instance, if the possibility of subgroups within samples is real, then this parameter can be used to refine the output by setting it to the expected size of the group. Note that the average values displayed in the results files are computed only on experiments where all thresholds could hold (referred to as “valid”). Hence an event with few valid experiments but high ratio differences will be placed above an event with many valid experiments but lower ratio differences.


Optional parameters - others
----------------------------

-h, --help: Display the help menu.

-j, --jxnonly: Flag to trigger the junction-only mode. If present, exon will not be considered in the analysis. This allows bypassing the execution of coverageBed which can be time-consuming. Note however, that without exons, retained introns cannot be identified. If set, the experimental design file does not need to contain any BAM or BED files. Three columns are sufficient: experiment, condition, and junction file.

-a, --annotation: GTF genome annotation file. The file is used to associate exons and genes to junctions. If no annotation is provided then the junction-only mode is activated. By default jSplice expects an Ensembl formatted GTF file but any other file can be used. One just has to specify the keyword for exon lines with the “-t” parameter.

-t, --type: Keyword to identify exon lines in the GTF annotation. Only those lines are used to create the genome. Default is “exon”, the Ensembl identifier.

-b, --nobam: Flag to indicate that coverageBed (or equivalent) has already been performed. Hence BED-like files are expected instead of BAM files in the experimental design file. Note that the genome used in later computation will be the one computed through the junctions. Hence if none of the exons given in the BED-like files matches the junctions, no exon information will be used. The BED-like files must obey the following rules:
•	A gene region must have the gene name or ID in the name field (4th column)
•	An exon must have “gene_name|EXN” in the name field. (The tag “|EXN” is used to distinguish genes from exons.)
•	The read count must be in the seventh column. (Similar to coverageBed output.)
If –j is set, then this parameter is irrelevant

-n, --nbcores: Set the number of cores/threads to use for coverageBed. By default, one thread is assigned per BAM file. If –j or –b is set, then this parameter is irrelevant.

-s, --samestrand: Flag for strand-specific RNA-seq (same strandedness). The flag is passed onto coverageBed.  That is, coverageBed will only report hits in B that overlap A on the same strand. If –j or –b is set, then this parameter is irrelevant.

-S, --diffstrand: Flag for strand-specific RNA-seq (different strandedness). The flag is passed onto coverageBed.  That is, coverageBed will only report hits in B that overlap A on the opposite strand. If –j or –b is set, then this parameter is irrelevant.

-e, --exon: Flag to indicate that junction files are exon-based. (E.g. Tophat junction files). Default is intron-based.
