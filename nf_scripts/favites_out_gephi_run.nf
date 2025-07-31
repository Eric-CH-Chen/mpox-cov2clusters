#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nextflow for nextstrain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : n/a
    Website: n/a
    Slack  : n/a
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    example run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nextflow run ../github-source/test_scripts/WIP/favites_out_gephi_run.nf -entry c_to_nstrain_single \
	--input_config data/favites_input.json \
	--output_dir runs/simulation_output --run_flag example_flag \
	-c config/nextflow_config

*/

/**
---------------------------------------------------------------------------------
Parameter definition
---------------------------------------------------------------------------------
*/
// Parameters - Must have
params.input_config
params.output_dir
params.auspice_config = "config/auspice_config_test.json"	// path con auspice config

params.test = ""

// Parameters - Optional
def date = new Date().format("yyyy-mm-dd")
//params.meta_tsv = ""	// path to meta file
params.thread_num = 1
params.run_id = "${date}"		// not required for now; uses run_flag
params.tags = "simulation"				// for when run_id isn't enough
params.help = false
params.title_str = "example run $params.run_flag - $params.run_id"
params.nextstrain = ""			// commands for nextstrain pipeline masking_ex_runs.nf
params.version = false
params.run_flag = "favitest-nf"

//params.auspice_config = "config/auspice_config_test.json"	// path con auspice config, for running nextstrain
params.module_nextstrain_path = 'nf_scripts/masking_ex_runs.nf'


/**
---------------------------------------------------------------------------------
Process import
---------------------------------------------------------------------------------
*/
include {nextstrain_favites_run_module} from params.module_nextstrain_path

// auxiliary scripts bin dir
params.bin_dir = workflow.launchDir + "bin/"

/**
---------------------------------------------------------------------------------
Header Info
---------------------------------------------------------------------------------
*/
def header() {
	
	log.info """
	Output data:	${params.output_dir}results/
	Output tree:	${params.output_dir}auspice/
	"""
}
def helpMe(){
	header()	// some param information
	exit 0 		// stop the script

}

process run_favites {
	publishDir "$params.output_dir", mode: 'copy'
	conda '/Users/ericchen/miniforge3/envs/favites64'

	input:
	path fl_config_file

	output:
	path "fav_out_dir/transmission_network.tsv", emit: network
	path "fav_out_dir/sequences.fas", emit: sample_seqs
	path "fav_out_dir/ancestral_sequence.fas", emit: ancestral_seq
	
	path "fav_out_dir"

	script:
	"""
	favites_lite.py -c $fl_config_file -o fav_out_dir
	"""
}

process proc_network_file {
	publishDir "$params.output_dir", mode: 'copy'
	conda '/Users/ericchen/miniforge3/envs/genomics64'
	
	input:
	path network_f

	output:
	path network_gephi_edge_f, emit: gephi_edge_network
	path network_gephi_node_name, emit: gephi_node_network
	path network_nx_meta_f, emit: nextstrain_meta

	script:
	network_gephi_edge_f = "n-gephi_edge_f.tsv"
	network_gephi_node_name = "n-gephi_node_time_f.tsv"
	network_nx_meta_f = "n-nx_f.tsv"


	"""
	# Prepare for gephi
	(echo "source\ttarget\ttime"; perl -pe 'BEGIN {\$count = 1} s/None/"None" . \$count++/ge' $network_f) > $network_gephi_edge_f

	# Prepare for nextstrain
	$workflow.launchDir/bin/favites_to_nx.R $network_f $network_nx_meta_f $network_gephi_node_name
	"""
}

process rename_favites_ancestral_seq {
	input:
	path favites_ancestral_seq

	output:
	path renamed_avites_ancestral_seq

	script:
	// prepare for nextstrain_run
	renamed_avites_ancestral_seq = "renamed_ancestral_file.fasta"
	"""
	perl -pe 's/^>.*\$/>Ancestral_Sequence/' $favites_ancestral_seq > $renamed_avites_ancestral_seq
	"""
}
process rename_favites_sampled_sequence {
	input:
	path favites_sample_seq

	output:
	path renamed_favites_sample_seq

	script:
	renamed_favites_sample_seq = "samples.fasta"
	"""
	perl $workflow.launchDir/bin/perl_rename_fav_seq.pl $favites_sample_seq $renamed_favites_sample_seq
	"""
}

workflow {
	proc_network_file(Channel.fromPath(params.test)) 
}

// Workflow from favites config file to nextstrain runs (via the 'masking_ex_runs.nf')
workflow c_to_nstrain_single {
	// Main Channel
	config_file_ch = Channel.fromPath(params.input_config)

	run_favites(config_file_ch)
	
	network_ch = run_favites.out.network.collect()

	// Process favites files
	//
	rename_favites_ancestral_seq(run_favites.out.ancestral_seq)
	rename_favites_sampled_sequence(run_favites.out.sample_seqs)
	proc_network_file(network_ch)

	nextstrain_favites_run_module(proc_network_file.out.nextstrain_meta, rename_favites_ancestral_seq.out, rename_favites_sampled_sequence.out, 
		Channel.fromPath(params.auspice_config))
	
}

workflow nextstrain_fav_run {
	// Main thing is that this particular run does not include the fix_tree.py step
	take:
	sample_seq
	ref_seq
	meta_file


	main:
	params.inputs_fasta = sample_seq
	params.reference_fasta = ref_seq
	params.meta_tsv = meta_file

	//include {test_run_favites} from params.module_nextstrain_path

	test_run_favites()
} 
