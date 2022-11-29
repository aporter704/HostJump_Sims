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
	**workflow up to here**
	- prune tips with removed reaction irrespective of location and export tree (constant sampling); use newick tree
	- prune tips with removed reaction for location 1 (humans) and all location 0 before first human sample. thereafter prune location 0 with some probability; use newick tree
	- prune sequence data set accordingly
	- estimate ml trees in iqtree
	- fit clock 
	- **Revise steps below**
	- prune tips tips with removed reaction location 1?
	- sample humans[2] with some rate REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE.tree ---> **python script 2**
	  - 4.1 sample reservoir with constant rate U[0.01, 0.5] REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE_1pconst_RESERVOIRCONSTANTRATE.tree
		- 4.2 sample reservoir after 1st human U[0.01, 0.05] REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE_1pskyline_RESERVOIRAFTERHUMANRATE.tree
	- 5. simulate sequences for 4. and then subsample for 4.1 and 4.2 ----> **shell command + R**
	- 6. estimate trees 4.1 and 4.2 ----> **shell command**
	- 7. for all trees calculate: number of imports | tmrca of first human import | tmrca of tree ---> **Rscript invoking dendropy for master tree. Can be done posthoc** 

  Potential implementations:
  	- sample migration rate from animal[1] to human[2] from U[a, b] -> MIGRATE ----> **python script 1**
    - consider merging removed and sampled compartments, because we will do sampling posthoc.
    
    
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
    
