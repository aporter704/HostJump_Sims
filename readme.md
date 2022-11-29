# Species jump project and to do list:

### Also check to do list on /sandbox
- workflow example for sars-cov-2:
	- make xml and name it something like REFID_2host_m12_MIGRATE.xml
	- run xml above to generate nexus tree file and json trajectories
	- export tree with locations annotated in tip labels
	- IMPORTANT TODOs:
		- ~~collapse singles before exporting trees!~~
		- ~~add sampling times to tip labels~~
		- ~~A script to pull putative importation events based on monophyletic groups!~~ 
	- simulate sequence data...
	- prune tips with removed reaction irrespective of location and export tree (constant sampling); use newick tree
	- prune tips with removed reaction for location 1 (humans) and all location 0 before first human sample. thereafter prune location 0 with some probability; use newick tree
	- prune sequence data set accordingly
	- estimate ml trees in iqtree
	- fit clock 
	- 7. for all trees calculate: number of imports | tmrca of first human import (branching off) | age of first human sample (age of first tip) | tmrca of tree ---> THIS NEEDS TO BE EDITED IN find_monophyletic.py 
    
    
    Rates and genome length:
    
    MPXV: 197,417 bp, 1e-5
    SARS-CoV-2: 29,850 bp, 9e-4 
    Ebola: 19,000 bp , 1.2e-3 
    HeV: 18,000bp, 6.5e-4 
    
    
    IQ-tree 2.0.3 & LSD2
    
    Using just alignment and dates: iqtree -s ALIGNMENT_FILE --date DATE_FILE(tab delimited)  

    Using alignment, dates and ML tree: iqtree -s ALIGNMENT_FILE --date DATE_FILE(tab delimited) -te TREE_FILE 
    
    
    --date-options -r a  -l
    Calculating CIs: --date-ci 100(resampling branch lengths 100 times) 
    Estimating root: -r a
    Extracting dates from taxon name after final "|" delimiter:  --date TAXNAME
    
