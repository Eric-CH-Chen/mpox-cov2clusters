#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nextflow for validating new params of cov2cluster for mpoxv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : n/a

----------------------------------------------------------------------------------------
Notes:
	Be mindful on which tree, meta-data, and reference-clust are used.
	cov2cluster depends on first two, only samples that exists in both are considered.
	 and the evaluation of cluster fit depends on the reference-clust (Shannon's low res).
	
	The text for outlier in cluster files are currently assume to be "-1" for outputs
	from cov2cluster, and "Outlier" from Shannon's manual overlay annotation. This will
	need to be adjusted in the cleanClustRef.R for newer input data set.
*/

/**
---------------------------------------------------------------------------------
Parameter definition
---------------------------------------------------------------------------------
*/
// Parameters - Must have
params.input_tree		// Newick tree
params.input_meta		// meta date of StrainID and Date
params.interval_file	// file of different combinations to test
params.steps = false	// to do automatic steps or not 
params.output_dir

// Parameters - Optional
def date = new Date().format("yyyy-mm-dd")
//params.benchRef = "data/Shannon_data/mpox_genomic_clusters_overlay_10May24.txt" // old data
params.benchRef = "data/Shannon_data/RawData/July_mpoxv/mpox_genomic_clusters_overlay_24July24.txt"
params.thread_num = 2
params.beta_str = ""				// specific beta string, only used for testing
params.cov2clust_param_default = ""
params.run_id = "${date}"			// not required for now; uses run_flag
params.tags = "run"					// for when run_id isn't enough
//params.title_str = "Run $params.run_flag - $params.run_id"
params.version = false
params.returnProb = "false"  // Passed into cov2clust to ask if it should return probability calculations

/**
---------------------------------------------------------------------------------
Example runs
---------------------------------------------------------------------------------

Simulation Run
nextflow run ../github-source/test_scripts/WIP/param-validation.nf -entry run_multiple  \
                --input_tree runs/simulation_output/simulation_raw.nwk \
                --input_meta runs/simulation_output/n-nx_f.csv \
                --output_dir runs/simulation_output/params/ \
				--benchRef runs/simulation_output/modularity_node_r6c2.csv \
                --interval_file data/exampleIntervalData.tsv

**/
// Process
process cov2cluster_step {
	//publishDir "$params.output_dir/tkout", mode: 'copy', saveAs: {filename -> "$runTag _step.txt"}
	publishDir "$params.output_dir/tkout", mode: 'copy', saveAs: {filename -> "${runTag}_${out_fstring}_${out_label_string}_step.txt"}
	//publishDir "$params.output_dir/tkout", mode: 'copy'
	tag "$runTag"
	input:
	path nwk_tree
	path dates_csv
	tuple val(runTag), val(co_int), val(co_branchD), val(co_dateD)
	val out_prob
	val out_fstring
	val out_label_string
	
	output:
	tuple path("*.txt"), val(runTag), emit: step

	script:
	beta_str = "c(" + co_int + ", " + co_branchD + ", " + co_dateD + ")"
	returnProb = params.returnProb
	
	// assign run label
	runLabel = ""
	if (out_label_string == "") {
		runLabel = runTag
	} else {
		runLabel = "BC"
	}
	"""
	$workflow.launchDir/bin/runCov2Cluster.R $nwk_tree $dates_csv '$beta_str' $out_prob $out_fstring $runLabel $returnProb $runTag
	"""
}

process cov2cluster_single {
	publishDir "$params.output_dir", mode: 'copy'

	input:
	path nwk_tree
	path dates_csv
	val beta_str
	val out_fstring
	val out_label_string
	val run_tag
	
	output:
	tuple path("*.txt"), val(run_tag), emit:  single

	script:
	def out_prob = 0.8
	returnProb = params.returnProb
	"""
	$workflow.launchDir/bin/runCov2Cluster.R $nwk_tree $dates_csv $beta_str $out_prob $out_fstring $out_label_string $returnProb
	"""
}

process clean_minus_1 {
	publishDir "$params.output_dir", mode: 'copy'

	// should only be used to remove '-1' for visualization
	input:
	path cluster_file
	
	output:
	path "*.tsv", emit: clean_cluster
	
	script:
	def clean_cluster_file = cluster_file.getSimpleName() + "_clean.tsv"
	"""
	#!/usr/bin/env -S Rscript 

	ori_cluster <- read.table("$cluster_file", header=TRUE, sep='\t', na.strings= '-1')

	clean_cluster <- ori_cluster[complete.cases(ori_cluster),]

	write.table(clean_cluster, file="$clean_cluster_file", sep = "\t", quote=FALSE, row.names=FALSE)

	"""
}

process cleanClust {
	publishDir "$params.output_dir/intermediate", mode: 'copy', pattern: "*.tsv"

	tag "$run_tag"

	/* 
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Note:
	This passes the modified cluster file, ie. Outlier -> Outlier1, to the next step,
	but saves only the unmodified cluster file to the publishDir
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	*/ 
	input:
	tuple path(run_clust_f), val(run_tag)
	path(ref_clust_f)

	output:
	tuple val(run_tag), path("*_clust.tsv"), path("*_ref.tsv"), emit: oriCleanCluster
	tuple val(run_tag), path("*_clust_rename.tsv"), path("*_ref_rename.tsv"), emit: cleanCluster
	

	script:
	"""
	$workflow.launchDir/bin/cleanClustRef.R $run_clust_f $ref_clust_f $run_tag
	"""

}

process calculateClusterFit {
	//publishDir "$params.output_dir/intermediate", mode: 'copy'
	// Assumes the cluster file and the reference cluster file has been cleaned and reordered!!
	// Input format (cluster_clean_file):
	//	col1	clust
	//	strain1	AUG18.NA.zzh
	// Input format (ref_clean_file):
	//	col1	clust
	//	strain1	clade1
	//
	// Output format, no header:
	//	clust_flag	sj_12	sj_21	adjRand_distance

	tag "$run_flag"
	input:
	tuple val(run_flag), path(cluster_clean_file), path(ref_clean_clust)

	output:
	//tuple val(run_flag), path("${clust_flag_f}"), emit: clustFitFlag // might not need this
	path("${clust_flag_f}"), emit: clustFit

	shell:
	clust_flag_f = run_flag +  ".csv"
	"""
	#!/usr/bin/env -S Rscript 
	library(mclust)
	library(igraph)

	r_flag <- "!{run_flag}" # declare run_flag in script

    # readClust
	clust_df <- read.table('!{cluster_clean_file}', sep="\t", header=TRUE, quote="")
	ref_df <- read.table('!{ref_clean_clust}', sep="\t", header=TRUE, quote="")
    
	# Get only the clustering clumn (#2)
	membership_1 <- clust_df[,2]
	membership_2 <- ref_df[,2]

	# adjRandIndex
	rIndx <- adjustedRandIndex(membership_1, membership_2)
	
	# split_join distanct
	# 	note: membership needs to be factors or integers
	sj_distance <- split_join_distance(as.factor(membership_1), as.factor(membership_2))

	# number of unclustered in clust
	#unclust_count <- length(which(membership_1 == "-1"))
	unclust_count <- length(which(grepl("^Outlier", membership_1) == TRUE))
	
	df_check <- data.frame("tag" = r_flag, "sj"=t(sj_distance), "sj_sum"= sum(sj_distance), "adjRand"=rIndx, "unclust_num" = unclust_count)
	
	write.table(df_check, file="!{clust_flag_f}", quote=FALSE, col.names=TRUE, row.names=FALSE, sep=",")
	"""
}

process collectAndPrepResultTable {
	publishDir "$params.output_dir/", mode: 'copy'
	input:
	path("*.csv")  		// list of multiple conditions
	path interval_file 	// for final file name

	output:
	path "*.tsv"

	script:
	"""
	$workflow.launchDir/bin/collectPrepResult.R $interval_file

	"""

}

// Workflows

workflow run_specific {
	tree = Channel.fromPath(params.input_tree).collect()
	meta = Channel.fromPath(params.input_meta).collect()
	ref = Channel.fromPath(params.benchRef).collect()
	def tagg = "Test"

	def outF = "mpox_" + params.run_id 
	cov2cluster_single(tree, meta, '"c(5,-3,-0.15,-0.2)"', outF, "single", tagg)

	cleanClust(cov2cluster_single.out.single, ref)
	calculateClusterFit(cleanClust.out.cleanCluster)
}

workflow run_multiple {
	// Raw input
	tree = Channel.fromPath(params.input_tree).collect()
	meta = Channel.fromPath(params.input_meta).collect()
	ref = Channel.fromPath(params.benchRef).collect()

	// Process interval file
	Channel.fromPath(params.interval_file)
		.splitCsv(header: true, sep: '\t')
		.map{ row -> tuple(row.RunTag, row.Int, row.BranchP, row.DateP) }
		.set{ interval_file_ch }

	// Run cov2cluster and collect outputs
	cov2cluster_step(tree, meta, interval_file_ch, 0.8, "mpox_Ex", "")

	// Clean outputs
	cleanClust(cov2cluster_step.out.step, ref)

	// Calculate clustering stats
	calculateClusterFit(cleanClust.out.cleanCluster)
	collectAndPrepResultTable(calculateClusterFit.out.clustFit.collect(), Channel.fromPath(params.interval_file))


}

