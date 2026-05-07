1. Do I need an interactive node? Is slurm the default?

Yes, MODE="slurm" is already the default. In slurm mode, Snakemake is just a lightweight orchestrator — it reads the Snakefile, builds the dependency graph, and calls sbatch to submit each rule as its own job. It does no computation itself, so you can run it from the login node. The key caveat is that if your SSH session drops, the Snakemake process dies (the already-submitted Slurm jobs keep running, but Snakemake won't submit the next ones). The fix is to wrap it in a screen session before launching:


screen -S dms_pipeline
cd /scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio
bash shell_scripts/run_dms_pipeline.bash
# Ctrl+A, D to detach; screen -r dms_pipeline to reattach
2. How does the pipeline know the inputs? Where does output go?

Everything is in config.yaml at the top of DMS_PacBio/:


fastq:    "...SMO.fastq.gz"     # input reads
genbank:  "...targets.gb"       # gene definitions
yaml_config: "...targets.yml"   # alignparse parsing parameters
outdir:   "...pipeline_out"     # ALL output lands here
For a new dataset in the future, copy and edit the config, then pass it explicitly:


snakemake -s Snakefile_library_qc --configfile config_newrun.yaml --executor slurm --jobs 50 ...
Or change CONFIG="config_newrun.yaml" at the top of the run script.

3. Are the scripts in stable locations?

The Python scripts you're relying on now live in /scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/python_scripts/, which is scratch and can be wiped. The originals are in /lts/rmlab/... (long-term storage). I'd recommend mirroring the working copy to lts once the pipeline is validated:


cp -r /scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/python_scripts/ \
      /lts/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/python_scripts_v2/
The Python environment (~/alignparse_env) and minimap2 (~/minimap2/) live in your home directory — those are stable and won't be wiped.

4. Naming — yes, make them specific

The Snakefile is now Snakefile_library_qc. The convention -s Snakefile_<name> lets you run any one of them from the same directory. When you build the screening pipeline next week, create Snakefile_screen_analysis alongside it. Each pipeline gets its own run wrapper (e.g., run_screen_analysis.bash) pointing to its specific Snakefile and config.