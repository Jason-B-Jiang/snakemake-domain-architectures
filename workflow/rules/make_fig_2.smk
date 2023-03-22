rule classify_lost_c_term_residues:
	input:
		'../../results/OrthoFinder/Results_OrthoFinder/Orthogroup_Sequences',
		'../../results/single_copy_orthogroups.csv',
		'../../data/sgd_to_uniprot_names.rds',
		'../../data/yeast_premature_stop_codons/supp_11.xls'
	output:
		'results/fig_2/lost_residue_distributions.svg',
		'results/fig_2/percent_dispensible_lost.svg',
		'results/fig_2/percent_losing_essential_residue.svg'
	shell:
		"Rscript scripts/classify_lost_c_terminal_residues.R {input} {output}"
