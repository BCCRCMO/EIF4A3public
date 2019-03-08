# EIF4A3public
Public repository for EIF4A3 paper with Takeda compounds

Publication: Pharmacological systems analysis defines EIF4A3 functions
in cell-cycle and RNA stress granule formation

Reference: Mazloomian et al (2019) Comms Biol


# Contents
Scripts that were used for the generation of the figures / analyses of the data

Figure 1.
	main/code/Gene_expression/cluster_genes.R
	main/code/Gene_expression/normalize.R

Figure 2.
	main/code/Miso/count_plots/make_plots.R
	main/code/Miso/count_plots/make_diff.py
	main/code/Vasttools/make_plots.R

Figure 3.
	main/code/Miso/clustering/cluster_ratios.R

Figure 4.
	main/code/Motifs/check_SE_distance.py
	main/code/Motifs/count_motifs.py
	main/code/Motifs/merged_statistical_test.py
	
Figure 5.
	main/code/Gene_expression/calculate_enrichment_for_sets.py
	
Figure 6.
	main/code/FACS/Caspase-time-dose.R
	
Figure 7.
	main/code/Stress_Granules/StressGranules.R
	
	
# Dependencies
In addition to the scripts posted, the following packages and software were used. 
These are not included in the repo, but available from the referenced sources.

Figure 1.
a - Adobe illustrator
b - ggplot package of R
citation info: https://cran.r-project.org/web/packages/ggplot2/citation.html
c - UpSetR package of R:
Alexander Lex, Nils Gehlenborg, Hendrik Strobelt, Romain Vuillemot, Hanspeter Pfister,
UpSet: Visualization of Intersecting Sets,
IEEE Transactions on Visualization and Computer Graphics (InfoVis '14), vol. 20, no. 12, pp. 1983–1992, 2014.
doi:10.1109/TVCG.2014.2346248

Figure 2.
a, b - ggplot package of R
c - matplotlip-venn package of python:
https://pypi.org/project/matplotlib-venn/

Figure 3.
a, b - ggplot package of R
c - googleVis package of R
https://cran.r-project.org/web/packages/googleVis/index.html

Figure 4.
a, b - seaborn package of python
https://seaborn.pydata.org/index.html
c - ggplot package of R

Figure 5.
a - EnrichmentMap (cited in the caption)
b, c - ggplot package of R

Figure 6.
a - ModFit software (Verity, Topsham, ME)
