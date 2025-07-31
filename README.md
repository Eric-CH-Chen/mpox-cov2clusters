# mpox-cov2clusters
In support of the manuscript submission.

# R libraries
  cov2clusters
  - ape
  - reshape2
  - rjson
  - stringi
  - data.table
  cleanClustRef
  - data.table
  - tools
  collectPrepResult
  - data.table
  - tools
  favites_to_nx
  - data.table
  - tools
  calculateClusterFit.nf
  - mclust
  - igraph

# General work flow
Step 1:
Use run favites_out_gephi_run.nf with a favites_lite.json file as input
ie. 
```
nextflow run ../github-source/test_scripts/WIP/favites_out_gephi_run.nf -entry c_to_nstrain_single \
	--input_config data/favites_input.json \
	--output_dir runs/simulation_output --run_flag example_flag \
	-c config/nextflow_config
```

