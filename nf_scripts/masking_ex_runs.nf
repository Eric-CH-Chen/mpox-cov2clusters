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

/**
---------------------------------------------------------------------------------
Parameter definition
---------------------------------------------------------------------------------
*/
// Parameters - Must have
params.reference_fasta
params.inputs_fasta
params.run_flag
params.meta_tsv	// path to meta file
params.output_dir
params.reference_genbank
params.auspice_config = "config/auspice_config_test.json"	// path con auspice config

// Parameters - Optional
def date = new Date().format("yyyy-mm-dd")
params.thread_num = 8
params.run_id = "${date}"		// not required for now; uses run_flag
params.tags = "run"	// for when run_id isn't enough
params.columns_header = "country datayear"	// example defaults
params.help = false
params.title_str = "example run $params.run_flag - $params.run_id"
params.version = false
params.aln_file = ""			// for when alignment already exists
params.ignore_sites	= false		// for when there are sites to ignore
params.ignore_beginning = false	// for when masking positions from the 5' end
params.ignore_end = false		// for when masking positions from the 3' end

params.align_fasta = 'null'			// aligned fasta, main reason being that it is a large alignment
//Parameters 

params.root = "NC_063383"		// used in fix_tree script
params.refine_root = "oldest"	// used in augur refine_tree; oldest is the default

// default using hmpox1 from nextstrain (commit: 27fdc7b, retrieved 2024-05-07)
params.clock_rate = 0			// fixed clock rate
params.clock_std_dev = 0		// standard deviation of fixed clock_rate estimate

// auxiliary scripts bin dir
params.bin_dir = workflow.launchDir

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

def version() {
  log.info """
  APOBEC3 analysis test run version: ${workflow.manifest.version}
  Nextflow version:	${nextflow.version}
  Command line:   ${workflow.commandLine}
  """
}
//version upon request
//if (params.version) {
//  version()
//  exit 0
//}
//displays help upon request
if (params.help) {
  helpMe()
  exit 0 //stop running
}
/**
---------------------------------------------------------------------------------
process definition
---------------------------------------------------------------------------------
*/
process nextstrain_mask {
	publishDir "$params.output_dir", mode: 'copy', saveAs: {filename -> "${params.run_flag}_aln_mask.fasta"}

	input:
	path align_fasta
	path exclude_sites
	val from_beginning
	val from_end

	output:
	path "mask_align.fasta", emit: mask_align_tree
	path "mask_align.fasta", emit: mask_align_fixtree
	path "mask_align.fasta", emit: mask_align_refine
	path "mask_align.fasta", emit: mask_align_ancestral
	

	script:
	def exclude_site_option = exclude_sites ? "--mask $exclude_sites" : ''
	def exclude_from_beinning_option = from_beginning ? "--mask-from-beginning $from_beginning" : ''
	def exclude_from_end_option = from_end ? "--mask-from-end $from_end" : ''
	"""
	augur mask --sequences $align_fasta --output mask_align.fasta \
		$exclude_site_option \
		$exclude_from_beinning_option \
		$exclude_from_end_option
	"""
}
process nextstrain_align {
	publishDir "$params.output_dir", mode: 'copy', saveAs: {filename -> "${params.run_flag}_align.fasta"}

	input:
	path input_fasta
	path reference_fasta

	output:
	path "align.fasta", emit: alignment_tree
	path "align.fasta", emit: alignment_fixtree
	path "align.fasta", emit: alignment_refine
	path "align.fasta", emit: alignment_ancestral

// example command to emulate
// augur align --nthreads -1 --sequences data/masking_example_set02/combined_ex_data.fasta --reference-sequence data/reference_genome_data/mpox_NextStrain_reference.fasta --output runs/example_mpox2/combined_ex_data_align.fasta  
	script:
	"""
	augur align --nthreads ${params.thread_num} --sequences $input_fasta \
		--reference-sequence $reference_fasta \
		--output align.fasta
	"""
}

process nextstrain_nextalign {
	publishDir "$params.output_dir", mode: 'copy', saveAs: {filename -> "${params.run_flag}_align.fasta"}
	// uses nextalign instead of augur align
	input:
	path input_fasta
	path reference_fasta

	output:
	path "align.fasta", emit: alignment_tree
	path "align.fasta", emit: alignment_fixtree
	path "align.fasta", emit: alignment_refine
	path "align.fasta", emit: alignment_ancestral

// Much faster, but a little different from augur align
	script:
	"""
	nextalign run -o align.fasta --include-reference --input-ref $reference_fasta $input_fasta

	"""
}

process nextstrain_nextclade_align {
	publishDir "$params.output_dir", mode: 'copy', saveAs: {filename -> "${params.run_flag}_align.fasta"}
	// uses nextalign instead of augur align
	input:
	path input_fasta
	path reference_fasta

	output:
	path "align.fasta", emit: alignment_tree
	path "align.fasta", emit: alignment_fixtree
	path "align.fasta", emit: alignment_refine
	path "align.fasta", emit: alignment_ancestral

// Much faster, but a little different from augur align
	script:
	"""
	nextclade3 run -o align.fasta  -r $reference_fasta --include-reference true $input_fasta

	"""
}

process nextstrain_tree {
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.nwk", saveAs: {filename -> "${params.run_flag}_raw.nwk"}
	input:
	path align_fasta
	path ignore_sites
	output:
	path "raw_tree.nwk", emit: raw_tree

// example command to emulate
// augur tree --alignment runs/example_mpox2/combined_ex_data_align.fasta --nthreads 8 \
//	--output runs/example_mpox2/combined_ex_data_rawtree.nwk

	script:
	def exclude_site_option = ignore_sites ? "--exclude-sites $ignore_sites" : ''
	"""
	augur tree --alignment $align_fasta --nthreads ${params.thread_num} \
		--output raw_tree.nwk $exclude_site_option
	"""
}

process reroot_tree {
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.nwk", saveAs: {filename -> "${params.run_flag}_reroot.nwk"}
	
	input:
	path raw_tree

	output:
	path "rerooted_tree.nwk", emit: reroot_tree

	script:
	"""
	gotree -h 
	"""

}
process nextstrain_refine {
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.nwk", saveAs: {filename -> "${params.run_flag}_refine.nwk"}
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.json", saveAs: {filename -> "${params.run_flag}_branch_length.json"}
	
	input:
	path input_raw_tree
	path align_fasta
	path meta_file

	output:
	// multiple duplicated channels for use in each downstream process
	path "refine_tree.nwk", emit: tree_export
	path "refine_tree.nwk", emit: tree_traits
	path "refine_tree.nwk", emit: tree_ancestral
	path "refine_tree.nwk", emit: tree_translate
	path "branch_length.json", emit: branch_length

// example command to emulate
// augur refine --tree runs/example_mpox2/combined_ex_data_rawtree.nwk \
//	--alignment runs/example_mpox2/combined_ex_data_align.fasta \
//	--metadata data/masking_example_set02/combined_ex_data_ini_meta.tsv \
//	--output-tree runs/example_mpox2/combined_ex_data_tree.nwk \
//	--output-node-data runs/example_mpox2/results/combined_ex_data_branch_length.json --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 4
	script:
	def clock_rate_param
	def clock_std_param

	if (params.clock_rate == 0 & params.clock_std_dev == 0) {
		// clock-rate and clock-std-dev not changed
	} else {

	}
	"""
	augur refine --tree $input_raw_tree \
		--alignment $align_fasta \
		--metadata $meta_file \
		--output-tree refine_tree.nwk \
		--output-node-data branch_length.json \
		--timetree --coalescent opt --date-confidence --date-inference joint \
		--keep-polytomies  \
		--divergence-units mutations \
		--precision 3 \
		--use-fft \
		--root $params.refine_root \
		
	"""
	// Changes:
	// 1. greedy-resolve is how things are currently done,removed per #4
	// 2. 'clock-filter-iqd 4' option, present in covid and zika tutorial, are removed here
	//	as we want to retain everything. All sequences are presumed good or interested.
	// 3. specify 'precision 3', which is the highest option,
	// 4. --keep-polytomies is used in both mpox and ncov, added here
	// 5. date-inference = joint in ncov, but marginal in mpox
	// 6. clock rate is now ignore or predefined. Eventually this should be redefined to be fed in 
//#--clock-rate 5.7e-5 --clock-std-dev 2e-5 
	// 7. clock rate reinstituted and seed removed
	// 8. -clock-rate 5.7e-5 --clock-std-dev 2e-5 removed for simulation


}

process nextstrain_traits {
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.json", saveAs: {filename -> "${params.run_flag}_traits.json"}

	input:
	path refine_nwk_tree
	path meta_file
	
	output:
	path "traits.json", emit: traits
	path "traits.json", emit: traits_export

// example command to emulate
// augur traits --tree runs/example_mpox2/combined_ex_data_tree.nwk \
//	--metadata data/masking_example_set02/combined_ex_data_fd_meta.tsv \
//	--output-node-data runs/example_mpox2/results/combined_ex_data_traits.json --columns country datayear --confidence
	// country may be taken out, as reconstruction of that isn't really needed.
	script:
	"""
	augur traits --tree $refine_nwk_tree \
		--metadata $meta_file \
		--output-node-data traits.json --columns ${params.columns_header} --confidence \
		--sampling-bias-correction 3
	"""
// --columns <country> <datayear>, this is based on meta data file, and correspond to the specific header
// will need to change later

}

process nextstrain_ancestral {
	// keeping this as the default
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.json", saveAs: {filename -> "${params.run_flag}_nt_muts.json"}

	input:
	path refine_nwk_tree
	path align_fasta
	
	output:
	path "nt_muts.json", emit: ancestral_export
	path "nt_muts.json", emit: ancestral_translate

// example command to emulate
// augur ancestral --tree runs/example_mpox2/combined_ex_data_tree.nwk \
//	--alignment runs/example_mpox2/combined_ex_data_align.fasta \
//	--output-node-data runs/example_mpox2/results/combined_ex_data_nt_muts.json --inference joint
	script:
	"""
	augur ancestral --tree $refine_nwk_tree \
		--alignment $align_fasta \
		--output-node-data nt_muts.json \
		--inference joint

	"""
// --inference joint, this is also from tutorial 
}

process nextstrain_ancestral_ignore_N {
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.json", saveAs: {filename -> "${params.run_flag}_nt_muts.json"}

	input:
	path refine_nwk_tree
	path align_fasta
	path ref_fasta
	
	output:
	path "nt_muts.json", emit: ancestral_export
	path "nt_muts.json", emit: ancestral_translate

	// Main difference being that it adds the option to ignore Ns and prevents snp calling on the overhangs
	// 	This should be the main one getting used.
	script:
	"""
	augur ancestral --tree $refine_nwk_tree \
		--alignment $align_fasta \
		--keep-ambiguous \
		--keep-overhangs \
		--output-node-data nt_muts.json \
		--root-sequence $ref_fasta \
		--inference joint

	"""
// --inference joint, this is also from tutorial 
}
process nextstrain_translate {
	tag "amino acid change estimate"
	publishDir "$params.output_dir", mode: 'copy', pattern: "*.json", saveAs: {filename -> "${params.run_flag}_aa_muts.json"}

	input:
	path refine_nwk_tree
	path reference_genbank
	path ancestral_json
	
	output:
	path "translate.json", emit: translate

// example command to emulate
// augur translate --tree runs/example_mpox2/combined_ex_data_tree.nwk \
//	--ancestral-sequences runs/example_mpox2/results/combined_ex_data_nt_muts.json \
//	--output-node-data runs/example_mpox2/results/combined_ex_data_aa_muts.json \
//	--reference-sequence data/reference_genome_data/mpox_NextStrain_reference.gb
	script:
	"""
	augur translate --tree $refine_nwk_tree \
		--ancestral-sequences $ancestral_json \
		--output-node-data translate.json \
		--reference-sequence $reference_genbank
	"""
}
process nextstrain_export {
	publishDir "$params.output_dir/auspice", mode: 'copy', pattern: "auspice.json", saveAs: {filename -> "${params.run_flag}_auspice.json"}
	publishDir "$params.output_dir/auspice", mode: 'copy', pattern: "auspice_root-sequence.json", saveAs: {filename -> "${params.run_flag}_auspice_root-sequence.json"}

	input:
	path refine_nwk_tree
	path node_jsons
	val title_str
	path meta_file
	path auspice_config

	output:
	path "auspice.json", emit: auspice
	path "auspice_root-sequence.json", optional: true

// example command to emulate
// ## NextStrain Export ##
// augur export v2 --tree runs/example_mpox2/combined_ex_data_tree.nwk \
//	--metadata data/masking_example_set02/combined_ex_data_fd_meta.tsv  \
//	--node-data _branch_length.json _traits.json _nt_muts.json a_aa_muts.json  \
//	--output runs/example_mpox2/auspice/combined_ex_data.json \
//	--title "mpox example set 2 raw" --maintainers Eric --auspice-config config/auspice_config_test.json
	script:
	"""
	## NextStrain Export ##
	augur export v2 --tree $refine_nwk_tree \
		--metadata $meta_file  \
		--node-data $node_jsons \
		--output "auspice.json" \
		--title "$title_str" --maintainers "Eric" --auspice-config $auspice_config \

	"""
	//		--include-root-sequence \
}
process collect_json {
	// To collect all json to feed into NextStrain export command
	//	Really is to make the script/commands more clear to read.
	input:
	path branch_length_json
	path traits_json
	path ancestral_json
	path translate_json

	output:
	path "*.json", includeInputs: true,  emit: jsons

	script:
	"""
	echo $branch_length_json $traits_json $translate_json
	"""
}

process py_genbank_convert {
	// not working yet; header issue
	// 
	input:
		path file_gbk
		val fasta_name

	output:
		path file_fasta
	
	def file_fasta = "$fasta_name"
	script:
	"""
	python -c "from Bio import SeqIO; SeqIO.convert($file_gbk, 'genbank', $file_fasta, 'fasta')"
	"""
}

process nx_fix_tree {
	//from nextstrain's mpox github

	input:
	path raw_tree
	path align_fasta
	val root_name

	output:
	path "tree_fixed.nwk", emit: fixed_tree

	script:
	def param_root = root_name ? "--root $root_name" : ""

	"""
	$workflow.launchDir/bin/fix_tree.py \
		--alignment $align_fasta \
		--input-tree $raw_tree \
		$param_root \
		--output "tree_fixed.nwk"

	"""
}
/**
---------------------------------------------------------------------------------
Workflow definition
---------------------------------------------------------------------------------
*/



workflow test_run_v3_nextclade {
	// WIP, be careful when running this workflow
	// Will be the main workflow, main difference being that "fix_tree added"
	// Set reference fasta and metadata tsv files as value channel
	ch_meta_tsv = Channel.fromPath(params.meta_tsv).collect()
	ch_ref_fasta = Channel.fromPath(params.reference_fasta).collect()

	nextstrain_nextclade_align(Channel.fromPath(params.inputs_fasta), ch_ref_fasta)
	
	// In case we need to mask, and this is done before tree-building, so these
	//	masked positions should be ignored in the subsequent stesps
	//	A more verbose way of declaring and checking if masking is requested
	def ig_sites = params.ignore_sites ? true : false
	def ig_beginning = params.ignore_beginning ? true : false
	def ig_end = params.ignore_end ? true : false

	def mask_param = [ig_sites, ig_beginning, ig_end]
	ch = Channel.of(mask_param)
	ch.view{"Masking: specific sites , 5' end, 3' end → $it"}

	// Initialize empty channel to save (masked) aligned fasta for tree, refine, and ancestral steps
	def ch_align_tree = []
	def ch_align_fixtree = []
	def ch_align_refine = []
	def ch_align_ancestral = []

	switch (mask_param) {
		case {it == [false, false, false]}:
			// no masking option declared 
			// go directly to tree making
			ch_align_tree = nextstrain_nextclade_align.out.alignment_tree
			ch_align_fixtree = nextstrain_nextclade_align.out.alignment_fixtree
			ch_align_refine = nextstrain_nextclade_align.out.alignment_refine
			ch_align_ancestral = nextstrain_nextclade_align.out.alignment_ancestral
			break;
		case {it == [true, false, false]}: 
			// only sites declared, mask the sites
			nextstrain_mask(nextstrain_nextclade_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), [], [])
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
		case {it == [true, true, true]}: 
			// mask both positions and the "ignore_beginning" positions from the star and "ignore_end" positions from the end
			nextstrain_mask(nextstrain_nextclade_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), Channel.of(params.ignore_beginning), Channel.of(params.ignore_end))
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
	}
	// make tree	
	nextstrain_tree(ch_align_tree, [])

	nx_fix_tree(nextstrain_tree.out.raw_tree, ch_align_fixtree, params.root)
	nextstrain_refine(nx_fix_tree.out.fixed_tree, ch_align_refine, ch_meta_tsv)
	// fix tree isn't working
	//nextstrain_refine(nextstrain_tree.out.raw_tree, ch_align_refine, ch_meta_tsv)

	// get annotation	
	nextstrain_traits(nextstrain_refine.out.tree_traits, ch_meta_tsv)
	nextstrain_ancestral_ignore_N(nextstrain_refine.out.tree_ancestral, ch_align_ancestral, ch_ref_fasta)
	nextstrain_translate(nextstrain_refine.out.tree_translate, Channel.fromPath(params.reference_genbank), nextstrain_ancestral_ignore_N.out.ancestral_translate)

	// collect jsons for export
	collect_json(nextstrain_refine.out.branch_length, 
					nextstrain_traits.out.traits_export, 
					nextstrain_ancestral_ignore_N.out.ancestral_export,
					nextstrain_translate.out.translate
				)

	// export to auspice for viewing
	nextstrain_export(nextstrain_refine.out.tree_export,
						collect_json.out.jsons.collect(), 
					params.title_str, ch_meta_tsv, 
					Channel.fromPath(params.auspice_config))
}
workflow test_run_favites{
	// WIP, be careful when running this workflow
	// This expects run using FAVITES output:
	//	- no alignment (prealigned)
	//  - assumes conversion of the input fasta
	ch_meta_tsv = Channel.fromPath(params.meta_tsv).collect()
	ch_ref_fasta = Channel.fromPath(params.reference_fasta).collect()

	nextstrain_align(Channel.fromPath(params.inputs_fasta), ch_ref_fasta)
	
	// In case we need to mask, and this is done before tree-building, so these
	//	masked positions should be ignored in the subsequent stesps
	//	A more verbose way of declaring and checking if masking is requested
	def ig_sites = params.ignore_sites ? true : false
	def ig_beginning = params.ignore_beginning ? true : false
	def ig_end = params.ignore_end ? true : false

	def mask_param = [ig_sites, ig_beginning, ig_end]
	ch = Channel.of(mask_param)
	ch.view{"Masking: specific sites , 5' end, 3' end → $it"}

	// Initialize empty channel to save (masked) aligned fasta for tree, refine, and ancestral steps
	def ch_align_tree = []
	def ch_align_fixtree = []
	def ch_align_refine = []
	def ch_align_ancestral = []

	switch (mask_param) {
		case {it == [false, false, false]}:
			// no masking option declared 
			// go directly to tree making
			ch_align_tree = nextstrain_align.out.alignment_tree
			ch_align_fixtree = nextstrain_align.out.alignment_fixtree
			ch_align_refine = nextstrain_align.out.alignment_refine
			ch_align_ancestral = nextstrain_align.out.alignment_ancestral
			break;
		case {it == [true, false, false]}: 
			// only sites declared, mask the sites
			nextstrain_mask(nextstrain_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), [], [])
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
		case {it == [true, true, true]}: 
			// mask both positions and the "ignore_beginning" positions from the star and "ignore_end" positions from the end
			nextstrain_mask(nextstrain_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), Channel.of(params.ignore_beginning), Channel.of(params.ignore_end))
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
	}
	// make tree	
	nextstrain_tree(ch_align_tree, [])

	//nx_fix_tree(nextstrain_tree.out.raw_tree, ch_align_fixtree, params.root)
	nextstrain_refine(nextstrain_tree.out.raw_tree, ch_align_refine, ch_meta_tsv)

	// get annotation	
	nextstrain_traits(nextstrain_refine.out.tree_traits, ch_meta_tsv)
	nextstrain_ancestral_ignore_N(nextstrain_refine.out.tree_ancestral, ch_align_ancestral, ch_ref_fasta)
	//nextstrain_translate(nextstrain_refine.out.tree_translate, Channel.fromPath(params.reference_genbank), nextstrain_ancestral_ignore_N.out.ancestral_translate)

	// collect jsons for export
	collect_json(nextstrain_refine.out.branch_length, 
					nextstrain_traits.out.traits_export, 
					nextstrain_ancestral_ignore_N.out.ancestral_export,
					[]
					//nextstrain_translate.out.translate
				)

	// export to auspice for viewing
	nextstrain_export(nextstrain_refine.out.tree_export,
						collect_json.out.jsons.collect(), 
					params.title_str, ch_meta_tsv, 
					Channel.fromPath(params.auspice_config))
}

workflow favites_run_params {
	// defines the params

	meta_file = params.meta_tsv
	reference_fasta = params.reference_fasta
	inputs_fasta = params.inputs_fasta
	auspice_config = params.auspice_config

	nextstrain_favites_run_module(meta_file, reference_fasta, inputs_fasta, auspice_config )


}

workflow nextstrain_favites_run_module {
	// defines the params
	take:
	meta_file_tsv
	reference_fasta
	inputs_fasta
	auspice_config

	main:
	// run scripts
	ch_meta_tsv = meta_file_tsv.collect()
	ch_ref_fasta = reference_fasta.collect()

	nextstrain_align(inputs_fasta, ch_ref_fasta)
	
	// In case we need to mask, and this is done before tree-building, so these
	//	masked positions should be ignored in the subsequent stesps
	//	A more verbose way of declaring and checking if masking is requested
	def ig_sites = params.ignore_sites ? true : false
	def ig_beginning = params.ignore_beginning ? true : false
	def ig_end = params.ignore_end ? true : false

	def mask_param = [ig_sites, ig_beginning, ig_end]
	ch = Channel.of(mask_param)
	ch.view{"Masking: specific sites , 5' end, 3' end → $it"}

	// Initialize empty channel to save (masked) aligned fasta for tree, refine, and ancestral steps
	def ch_align_tree = []
	def ch_align_fixtree = []
	def ch_align_refine = []
	def ch_align_ancestral = []

	switch (mask_param) {
		case {it == [false, false, false]}:
			// no masking option declared 
			// go directly to tree making
			ch_align_tree = nextstrain_align.out.alignment_tree
			ch_align_fixtree = nextstrain_align.out.alignment_fixtree
			ch_align_refine = nextstrain_align.out.alignment_refine
			ch_align_ancestral = nextstrain_align.out.alignment_ancestral
			break;
		case {it == [true, false, false]}: 
			// only sites declared, mask the sites
			nextstrain_mask(nextstrain_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), [], [])
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
		case {it == [true, true, true]}: 
			// mask both positions and the "ignore_beginning" positions from the star and "ignore_end" positions from the end
			nextstrain_mask(nextstrain_align.out.alignment_tree, Channel.fromPath(params.ignore_sites), Channel.of(params.ignore_beginning), Channel.of(params.ignore_end))
			ch_align_tree = nextstrain_mask.out.mask_align_tree
			ch_align_fixtree = nextstrain_mask.out.mask_align_fixtree
			ch_align_refine = nextstrain_mask.out.mask_align_refine
			ch_align_ancestral = nextstrain_mask.out.mask_align_ancestral
			break;
	}
	// make tree	
	nextstrain_tree(ch_align_tree, [])

	//nx_fix_tree(nextstrain_tree.out.raw_tree, ch_align_fixtree, params.root)
	nextstrain_refine(nextstrain_tree.out.raw_tree, ch_align_refine, ch_meta_tsv)

	// get annotation	
	nextstrain_traits(nextstrain_refine.out.tree_traits, ch_meta_tsv)
	nextstrain_ancestral_ignore_N(nextstrain_refine.out.tree_ancestral, ch_align_ancestral, ch_ref_fasta)
	//nextstrain_translate(nextstrain_refine.out.tree_translate, Channel.fromPath(params.reference_genbank), nextstrain_ancestral_ignore_N.out.ancestral_translate)

	// collect jsons for export
	collect_json(nextstrain_refine.out.branch_length, 
					nextstrain_traits.out.traits_export, 
					nextstrain_ancestral_ignore_N.out.ancestral_export,
					[]
					//nextstrain_translate.out.translate
				)

	// export to auspice for viewing
	nextstrain_export(nextstrain_refine.out.tree_export,
						collect_json.out.jsons.collect(), 
					params.title_str, ch_meta_tsv, 
					auspice_config)
}

workflow.onComplete {

    def msg = """\
        
        APOBEC3 Testing scripts
        
        Pipeline Execution Summary
        --------------------------------
        Completed at:   ${workflow.complete}
        Duration:   ${workflow.duration}
        Success:    ${workflow.success}
        Command line:   ${workflow.commandLine}
        Output Data directory:   ${params.output_dir}
        Nextflow version:   ${nextflow.version}
        Exit status:    ${workflow.exitStatus}
        """
        .stripIndent()

        log.info msg
}