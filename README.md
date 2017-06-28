# jsplice 1.0.2
jSplice is a fast method to detect differential alternative splicing events from RNA-seq
experiments. (see publication and/or manual for details).

jSplice is freely available under the GPL license at http://www.mhs.biol.ethz.ch/research/krek/jsplice

Requirements
------------
- Python 2.7 with numpy/scipy modules
- coverageBed from BEDtools package


Installation
------------
Nothing required

Execution
---------
Run 'python run.py -h' (or 'python -m jsplice.run' if the jsplice directory is in your PYTHONPATH) to access the help menu

What's new in version 1.0.2
---------------------------
- Bug fixes
- Added sanity checks on RPKM calculations
- 3-decimal values in HTML output
- Better handling of Ensembl and UCSC chromosome naming schemes
