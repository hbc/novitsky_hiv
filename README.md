## HIV analysis pipeline

Data, scripts and documentation for running HIV analysis:

- Quality control of human/HIV with bcbio runs
- Assembly and minority variant detection with existing pipelines from the
  [ICONIC](https://www.ucl.ac.uk/health-informatics/research/impact-in-research/iconic) project:
  [ICONIC-URL](https://github.com/ICONIC-UCL)

## Running

### Alignment and quality control with bcbio

Runs an initial alignment against HIV and human genome build 37, disambiguating
reads and assigning them to either. Then runs quality control metrics and
prepares an overall HTML report. There is a lot of documentation about the
project and usage: http://bcbio-nextgen.readthedocs.io/en/latest/

1. Add bcbio_nextgen.py to your PATH. Practically, we maintain a bcbio
   installation on Odyssey, include it with:

        export PATH=/n/regal/hsph_bioinfo/bcbio_nextgen/bin:$PATH

2. Create a bcbio configuration file. This happens in 3 steps:

   a. Make a CSV file describing the samples. There is an example CSV file in
      the config directory on Odyssey for the first set of 96 samples.
   b. Create a template YAML configuration file describing what to run. There is
      also an example file in the config directory.
   c. Run `bcbio_nextgen.py -w template` on the sample CSV, template and set of
      fastq files. There is a script in the config directory that does this.

3. Run analysis from configuration file, distributing on the cluster. We
   normally run from a `work` directory and there is an example script on
   Odyssey for how to submit. In non-Odyssey environments you can run in
   parallel on a single machine with multiple cores using:

        bcbio_nextgen.py ../config/your_config.yaml -n 8`

### UCL pipeline

1. Need a ready to run UCL installation with data and scripts. This is ready to
   go on Odyssey, see path documentation below.

2. There are some options to edit at the top of
   `ICONIC_UCL_pipeline/Scripts/command_line_HIV_pipeline.sh`. You should not
   need to change these for a ready to go installation like on Odyssey.

3. Create tab-delimited file of sampleID, forward read path, reverse read path.
   This should be called `samples.tsv` and can be in any directory.

4. Run job for each sample number in the file. There are example scripts on
   Odyssey that do this for multiple samples, but the basic command line is:

        bash ICONIC_UCL_pipeline/Scripts/command_line_HIV_pipeline.sh 1

   For a practical run, copy the `run.sh` script from an `iconic_ucl_pipeline`
   directory and edit to specify the samples you want run.

By default the pipeline uses 8 cores for each sample run. It ran relatively
quickly for 96 samples using 8 cores today on a single machine.

## Odyssey details

### Base directories

- Base: `/n/chb/projects/novitsky_hiv`
- UCL Pipeline code: `ICONIC_UCL_pipeline`
- UCL Pipeline runs: `FLOWCELL/iconic_ucl_pipeline/`

### UCL pipeline Example runs

- All samples: `FLOWCELL/iconic_ucl_pipeline/run_all`
   - samples.tsv -- list of all samples
   - run.sh -- script to submit job (`sbatch run.sh`)

- Single sample: `FC02539/iconic_ucl_pipeline/GEN00075776`
   - samples.tsv -- single sample list
   - run.sh -- script to submit job (`sbatch run.sh`)

You need to supply your own `samples.tsv` for each run, and can generally use
the `run.sh` script with only changes to the sample numbers to run.

### QC pipeline

Run in the base directory using bcbio (https://github.com/chapmanb/bcbio-nextgen)

- Configuration files: `FC02668/config/FC02668.yaml`
   - Script to create sample CSV file from fastqs: `prep_samples.py`
   - Script to create final configuration: `prep_config.sh`
- Run analysis: work
   - `cd work && sbatch ../submit_bcbio.sh`
   - final output QC report: `FC02668/final/2016-11-08_FC02539/multiqc/multiqc_report.html`

You can copy and use these three scripts but will have to adjust input files and output
naming to match the new flowcell.

## Installation

Unpack script code and Data (todo, put this into GitHub for automated download,
coordinate with Dan)

Install dependencies with [bioconda](https://bioconda.github.io/):

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p tools
    ./tools/bin/conda install -y -c conda-forge -c bioconda iva lastz trimmomatic blast parallel
