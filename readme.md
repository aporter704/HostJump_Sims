# Species jump project and to do list:

### Also check to do list on /sandbox
- workflow example for sars-cov-2:
	- make xml and name it something like REFID_2host_m12_MIGRATE.xml
	- run xml above to generate nexus tree file and json trajectories
	- sample humans[2] with some rate REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE.tree ---> **python script 2**
	  - 4.1 sample reservoir with constant rate U[0.01, 0.5] REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE_1pconst_RESERVOIRCONSTANTRATE.tree
		- 4.2 sample reservoir after 1st human U[0.01, 0.05] REFID_2host_m12_MIGRATE_2pHUMANSAMPRATE_1pskyline_RESERVOIRAFTERHUMANRATE.tree
	- 5. simulate sequences for 4. and then subsample for 4.1 and 4.2 ----> **shell command + R**
	- 6. estimate trees 4.1 and 4.2 ----> **shell command**
	- 7. for all trees calculate: number of imports | tmrca of first human import | tmrca of tree ---> **Rscript invoking dendropy for master tree. Can be done posthoc** 

  Potential implementations:
  	- sample migration rate from animal[1] to human[2] from U[a, b] -> MIGRATE ----> **python script 1**
    - consider merging removed and sampled compartments, because we will do sampling posthoc.
