=================================
Mouse Cortex + Hippocampus
=================================

RNA sequencing data of single cells isolated from >20 areas of mouse cortex and hippocamopus, including ACA, AI, AUD, CA, CLA, CLA;EPd, ENTl, ENTm, GU;VISC;AIp, HIP, MOp, MOs, ORB, PAR;POST;PRE, PL;ILA, PTLp, RSP, RSPv, SSp, SSs, SSs;GU, SSs;GU;VISC, SUB;ProS, TEa;PERI;ECT, VISal;VISl;VISli, VISam;VISpm, VISp, and VISpl;VISpor.  Abbreviations match the Allen Mouse Brain Atlas.

The data set includes 76,307 single cells.
The sequencing results were aligned to exons and introns in the GRCm38.p3 reference genome using the STAR algorithm, and aggregated intron and exon counts at the gene level were calculated.
For more details, please see the Documentation tab in the Cell Types web application.


Gene expression data matrices (transcriptome.zip)
    This file is a wrapper for several sparse matrices of gene expression values in HDF5 format, as well as some additional accessory pieces of information in other formats.
	/data/exon
		Contains the (row, column) matrix of read counts obtained for each (cell, gene) based on alignment to exons.
	/data/intron
		Contains the (row, column) matrix of read counts obtained for each (cell, gene) based on alignment to introns.
	/data/t_exon
		Contains the (row, column) matrix of read counts obtained for each (gene, cell) based on alignment to exons.
	/data/t_intron
		Contains the (row, column) matrix of read counts obtained for each (gene, cell) based on alignment to introns.
	/sample_names
	    Unique identifiers of the samples (cells).  This can be matched to the "sample_name" column in other files.
	/gene_names
	    Unique identifiers of the genes.  This can be matched to the "gene_symbol" column in other files.
	
	"tome" is a format developed at the Allen Institute with the primary goal to combine compact storage with reasonably fast random access of both genes and samples.
	The "scrattch.io" R library was developed for input, manipulation, and output of tome format (https://github.com/AllenInstitute/scrattch.io).
	Specific instructions for accessing the above data are as follows:
		1. Install the R programming language, as described here: https://cran.r-project.org/
		2. Install the "scrattch.io" library and its required dependencies, as described here: https://github.com/AllenInstitute/scrattch.io
		3. Navigate to the folder where you downloaded "transcrip.tome" and start R (or start R and then use "setwd" to navigate to the folder).
		4. Read in the relevant data matrix using these lines of code:
				library(scrattch.io)
				options(stringsAsFactors = FALSE)
				tome        <- "transcrip.tome"
				exons       <- read_tome_dgCMatrix(tome,"data/t_exon")    # (or data/exon)
				introns     <- read_tome_dgCMatrix(tome,"data/t_intron")  # (or data/intron)
				sample_name <- read_tome_sample_names(tome)  
				gene_name   <- read_tome_gene_names(tome)
				# Note that the sample and gene names are NOT stored within the exon and intron matrices and must be read separately. 
		5. IMPORTANT NOTE: all other matrices contained within "transrip.tome" are archival and SHOULD BE IGNORED.

		
Medians (medians.zip)
	A table of median expression values for each gene (rows) in each cluster (columns).  Medians are calculated by first normalizing gene expression as follows: norm_data = log2(CPM(exons+introns)), and then calculating the medians independently for each gene and each cluster.
	The first row lists the cluster name (cluster_label), which matches the cell type alias shown in the Transcriptomic Explorer.
	The first column lists the unique gene identifier (gene), which in most cases is the gene symbol.
	
	
Probability matrix of cluster membership (sample_cluster_probabilities.zip)
	A table indicating the probability that each cell (rows) is a member of each cell type (columns).
	Values are derived using a bootstrapping approach, where each cell is mapped to the taxonomy 100 times using different subsets of the data and the fraction of times each cell is mapped to each cluster is stored as the probability.  For more details on this method, please see Tasic et al 2018 (https://www.nature.com/articles/s41586-018-0654-5).
	The first row lists the cell type accession number (cell_type_accession_label) for each cluster.
	The first column lists the unique sample identifier (sample_name) for each cell. 


Cell metadata ("sample-annotations.zip")
* Each item of this table (except "sample_name") has three columns:
	[item]_label
		Name of the item (e.g., "V1C" would be an example under "brain_region_label")
	[item]_order
		Order that the item will be displayed on the Transcriptomics Explorer 
	[item]_color
		Color that the item will be displayed on the Transcriptomics Explorer 

* Items in the sample information table:
	sample_name
		Unique sample identifier
	cluster
		Cell type cluster name
	cell_type_accession
		Cell type accession ID (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	cell_type_alias
		Cell type alias (see https://portal.brain-map.org/explore/classes/nomenclature for details).  This is the same as "cluster".
	cell_type_alt_alias
		Cell type alternative alias, if any (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	cell_type_designation
		Cell type label (see https://portal.brain-map.org/explore/classes/nomenclature for details)
	class
		Broad cell class (for example, "GABAergic", "Non-neuronal", and "Glutamatergic")
	subclass
		Cell type subclass (for example, "SST", "L6 CT", and "Astrocyte")
	external_donor_name
		Unique identifier for each mouse donor
	donor_sex
		Biological sex of the donor
	cortical_layer
		Cortical layer targeted for sampling
	region
		Brain region targeted for sampling
	subregion
		Brain sub-region targeted for sampling (e.g., anterior vs. posterior), if any
	full_genotype
		Full genotype of the transgenic mouse donor
	facs_population_plan
		FACS gating criteria used to sort labeled cells
	injection_materials
		Specific virus injected into the mouse.  Blank values for this and subsequent columns indicate that no injection was performed.
	injection_method
		Method used for virus injection (Nanoject, Retro-Orbital)
	injection_roi
		Center of injection site. Abbreviations match the Allen Mouse Brain Atlas.
	propagation_type
		Type of viral propogation (retrograde, anterograde)
	
	
TSNE coordinates (2d-coordinates.zip)
t-Distributed Stochastic Neighbor Embedding (t-SNE) coordinates for each sample shown on the Transcriptomics Explorer.  t-SNE is a method for dimensionality reduction of gene expression that is  well suited for data visualization (as of 1 October 2019, a comprehensive t-SNE resource is available here: https://lvdmaaten.github.io/tsne/)
	sample_name
		Unique sample identifier
	tsne_1
		First t-SNE coordinate
	tsne_2
		Second t-SNE coordinate

		
Taxonomy of clusters (dendrogram.zip)
	Serialized cluster hierarchy with all node information embedded in json format.
	The dendrogram shown at the top of the Transcriptomics Explorer, including the underlying cluster order, is derived from this file.
	

Taxonomy metadata (taxonomy.zip)
	Tracking taxonomy meta-data is critical for reproducibility.  This file is a draft of taxonomy meta-data to be stored.  See the "Tracking taxonomies" section at https://portal.brain-map.org/explore/classes/nomenclature for details of each descriptor.
	

Gene information (**STORED ELSEWHERE**)
* To access this file, please use the following link: http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985
* Within that zip file, the gene information is located in "mouse_VISp_2018-06-14_genes-rows.csv".  All other files can be ignored.
	gene_symbol
		Gene symbol
	gene_id
		This is an Allen Institute gene ID that can be ignored
	chromosome
		Chromosome location of gene
	gene_entrez_id
		NCBI Entrez ID
	gene_name
		Gene name

		
Gene ".gtf" file (**STORED ELSEWHERE**)
* To access this file, please use the following link: http://celltypes.brain-map.org/api/v2/well_known_file_download/502999254
.gtf is a standard format for localizing various aspects of transcripts within a specific genome and information about this format is plentiful.
As of 1 October 2019, one active link describing this format is here: https://www.gencodegenes.org/pages/data_format.html