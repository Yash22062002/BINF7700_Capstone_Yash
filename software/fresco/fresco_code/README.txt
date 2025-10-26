				FRESCO (Finding Regions of Excess Synonymous Constraint)
				UNIX/Linux/Mac OS X (command-line)
				Last updated: May 27, 2014

Outline:
A. INTRODUCTION
B. RUNNING FRESCO
C. SAMPLE INPUT FILES
D. OUTPUT
E. FAQ
F. LINKS AND SUPPORT
G. CITATION
———————————————————————————————————————————————————————————————————
A. INTRODUCTION

FRESCO (Finding Regions of Excess Synonymous COnstraint) is a method designed to identify regions with excess synonymous constraint in short, deep alignments. Regions of excess synonymous constraint are common in the genomes of many types of organisms, and may indicate the presence of overlapping functional elements.

We provide code for running this method, implemented as a HYPHY batch script. You will first need to download and install the command-line version of the HYPHY package to run this script.
———————————————————————————————————————————————————————————————————
B. RUNNING FRESCO

1. Download and install the command line version of the HYPHY package (Pond et al., 2005). As of the time of writing this help file, you can find this at: http://hyphy.org/w/index.php/Download

2. Prepare the input files for the dataset that you are interested in. FRESCO requires two input files:
	(a) A sequence alignment in fasta format
	(b) A tree in Newick format
    NOTE that only alphanumeric characters and underscores are allowed in sequence names, all sequence names must be unique, and the names of sequences in the tree and the alignment must exactly match.

3. Edit runFRESCO.bf. This is a helper batch file that provides the input arguments to FRESCO. You will need to specify four things:
	(a) The file path to your sequence alignment
	(b) The file path to your phylogenetic tree
	(c) The size of the sliding window (in codons)
	(d) The file path to the FRESCO.bf batch file.
    NOTE that you will need to provide the either the complete path to your sequence alignment and tree, or the path relative to the FRESCO.bf file.

4. You can then run FRESCO using the following command:

/path/to/HYPHYMP runFRESCO.bf

(where you have replaced /path/to/ with the correct path to the HYPHY executable on your system).
———————————————————————————————————————————————————————————————————
C. SAMPLE INPUT FILES

As example input files, we have included the following:
	(1) HBV.P.10.1.fa: Codon alignment of the polymerase gene of ten isolates of Hepatitis B virus.
	(2)RAxML_bestTree.HBV.P.10.1.tree: Maximum-likelihood tree constructed using RAxML (using the GTRGAMMA nucleotide substation model and the above codon alignment).
———————————————————————————————————————————————————————————————————
D. OUTPUT

FRESCO will print to standard output a table containing information about excess synonymous constraint in the data. The output table has the following columns:
	(1) window_start: The codon position of the first codon in the sliding window (0-based indexing).
	(2) syn_rate_alternate_model: The synonymous rate under the alternate model (window-specific synonymous rate).
	(3) syn_rate_null_model: The synonymous rate under the null model (no window-specific synonymous rate; this column is always 1).
	(4) nonsynonymous_rate_alternate_model: The non-synonymous rate under the alternate model.
	(5) nonsynonymous_rate_null_model: The non-synonymous rate under the null model.
	(6) log_likelihood_alternate_model: The log-likelihood of the alternate model
	(7) log_likelihood_null_model: The log-likelihood of the null model
	(8) n_local_params_alternate_model: number of local parameters fit to each sliding window in the alternate model; this column is always 2.
	(9) n_local_params_null_model: number of local parameters fit to each sliding window in the null model; this column is always 1.
	(10) LRT: Result of the likelihood ratio test of the null and alternate models.
	(11) p_value: Gives the (uncorrected) p-value that the alternate model is a significantly better fit to the data within the given window than the null model.

By importing the output table into a data visualization environment (for example, R), you can easily plot the local synonymous and non-synonymous substitution rate across your sequence alignment.

    NOTE that the p_values in column 11 are uncorrected. You will want to correct for multiple hypothesis testing before reporting regions of significant excess synonymous constraint.



———————————————————————————————————————————————————————————————————
E. FAQ:

Q. I see an error regarding my sequence / tree format?

A. FRESCO is written using the HYPHY batch language (Pond et al., 2005), so names of the input sequences must be valid in the HYPHY language. (Only alphanumeric characters and underscores are allowed in sequence names and all sequence names must be unique). Also, the sequence names in the provided phylogenetic tree must perfectly match the sequence names in the alignment. Note also that all nodes in the provided Newick tree are required to have unique names (internal nodes do not need names, but naming internal nodes with non-unique names corresponding to bootstrap values may cause an error).

You may ALSO see this error if you have incorrectly specified the path to your sequence or tree file. Check that the path to your sequence and tree files are complete and correct.

Q. What do I need to provide as input to FRESCO?

A. You will first need to construct a codon alignment and a phylogenetic tree. We recommend first aligning the amino acid sequence of your proteins and then using the aligned proteins as a guide for a codon alignment. FRESCO reads alignments in FASTA format and trees in Newick format.

Q. How should I choose the window size?

A. Different types of constraint are best detected at different window sizes. Strong, short regions of constraint may be best pinpointed by smaller windows, while longer, weaker signals of excess synonymous constraint may be detected only at larger window sizes. We suggest trying a range of window sizes (for example, 1, 5, 10, 20 and 50-codon sliding windows).
———————————————————————————————————————————————————————————————————
F. SUPPORT

Please let us know if you find bugs or have suggestions for improvement: 

Rachel Sealfon
rsealfon@mit.edu
———————————————————————————————————————————————————————————————————
G. CITATION

Please cite the FRESCO paper when using or adapting this code.
