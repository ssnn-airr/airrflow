# nf-core/airrflow: Single-cell tutorial

This tutorial provides a step by step introduction on how to run nf-core/airrflow on single-cell BCR-seq data or single-cell TCR-seq data.

## Pre-requisites

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow and a container engine needed to run this pipeline. At the moment, nf-core/airrflow does NOT support using conda virtual environments for dependency management, only containers are supported. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) before running the workflow on actual data.

For the purpose of running this tutorial on your local machine, we recommend a docker installation. To install docker, follow the instructions [here](https://docs.docker.com/engine/install/). After docker installation on linux system, don't forget to check the [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/).

Alternatively, you can run this tutorial using the Gitpod platform which comes with Nextflow, nf-core and docker pre-installed. Please follow [these instructions](https://nf-co.re/docs/tutorials/gitpod/overview) on how to set it up.

## Testing the pipeline with built-in tests

Once you have set up your Nextflow and container (docker or singularity), test the airrflow pipeline with built-in test.

```bash
nextflow run nf-core/airrflow -r 4.2.0 -profile test,docker --outdir test_results
```

If the tests run through correctly, you should see the execution of airrflow processes. Finally, the following output will appear in your command line:

```bash
output:

-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 11:30:35
Duration    : 5m 50s
CPU hours   : 0.6
Succeeded   : 221
```

## Supported input formats

There are two supported input formats for nf-core/airrflow single-cell AIRRseq pipeline: assembled sequences in AIRR rearrangement format or departing from raw reads in fastq format sequenced in the 10x Genomics platform.

The [AIRR rearrangement format](https://docs.airr-community.org/en/latest/datarep/rearrangements.html) is a standard format to store BCR and TCR sequence data with relevant metadata fields. This format is supported as input and output by multiple tools specific for analyzing AIRR sequencing data. For example, when analyzing single-cell AIRR sequencing data with CellRanger versions >= 4.0 an AIRR rearrangement file will be provided as output, and this is the recommended input for running nf-core/airrflow. Note that it is also possible to start running the pipeline directly from raw sequencing reads, and in this case CellRanger will be run when launching nf-core/airrflow.

The AIRR rearrangement format is also the default one when analyzing publicly available data from specialized AIRRseq databases such as the AIRR Data Commons through the [iReceptor gateway](https://gateway.ireceptor.org/login).

In this tutorial we will showcase how to run nf-core/airrflow with both of the input formats.

![nf-core/airrflow overview](https://raw.githubusercontent.com/nf-core/airrflow/master/docs/images/airrflow_workflow_overview.png)

## Starting from AIRR rearrangement format

### Datasets

For this tutorial we will use subsampled PBMC single-cell BCR sequencing data from two subjects, before (d0) and after flu vaccination (d12).
The dataset is publicly available on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.11373740).
You don't need to download the dataset bacause the links to the samples are already provided in the samplesheet and Nextflow will get the data from the links automatically when running the pipeline.

### Preparing the samplesheet and configuration file

To run the pipeline, a tab-separated samplesheet that provides the path to the AIRR rearrangement files must be prepared.
The samplesheet collects experimental details that are important for the data analysis.
Details on the required columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#assembled-input-samplesheet-bulk-or-single-cell-sequencing).

The resource configuration file sets the compute infrastructure maximum available number of CPUs, RAM memory and running time. This will ensure that no pipeline process requests more resources than available in the compute infrastructure where the pipeline is running. The resource config should be provided with the `-c` option. In this example we set the maximum RAM memory to 16GB, we restrict the pipeline to use 8 CPUs and to run for a maximum of 24 hours.

```json title="resource.config"
process {
    resourceLimits = [ memory: 16.GB, time: 24.h, cpus: 8 ]
}
```

A prepared samplesheet for this tutorial can be found [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/assembled_samplesheet.tsv), and the configuration file is available [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/resource.config).
Download both files to the directory where you intend to run the airrflow pipeline.

> [!TIP]
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more cpus than available or run out of memory. Depending on the size of your dataset, it might be required to extend the running time. You can also remove the "time" parameter from the configuration file to allow for unlimited runtime.

> [!NOTE]
> When running nf-core/airrflow with your own data, provide the full path to your input files under the filename column.

### Running airrflow

With all the files ready, you can proceed to start the pipeline run:

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode assembled \
--input assembled_samplesheet.tsv \
--outdir sc_from_assembled_results  \
-c resource.config \
-resume
```

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/airrflow_sc_from_assembled.sh).
With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_sc_from_assembled.sh
```

> [!TIP]
> When launching a Nextflow pipeline with the `-resume` option, any processes that have already been run with the exact same code, settings and inputs will be cached and the pipeline will resume from the last step that changed or failed with an error. The benefit of using "resume" is to avoid duplicating previous work and save time when re-running a pipeline.
> We include "resume" in our Nextflow command as a precaution in case anything goes wrong during execution. After fixing the issue, you can relaunch the pipeline with the same command, it will resume running from the point of failure, significantly reducing runtime and resource usage.

After launching the pipeline the following will be printed to the console output, followed by some Nextflow parameters and executions of Airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

WARN: It appears you have never run this project before -- Option `-resume` is ignored
Launching `https://github.com/nf-core/airrflow` [boring_heyrovsky] DSL2 - revision: d91dd840f4 [4.2.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 4.2.0
------------------------------------------------------

```

Once the pipeline has finished successfully, the following message will appear:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 13:06:05
Duration    : 2m 47s
CPU hours   : 0.4
Succeeded   : 44
```

## Starting from raw reads in fastq format

### Datasets

For this tutorial we will use subsampled blood single-cell TCR sequencing data of one subject generated from the 10x Genomic platform. The links to the fastq files are in the samplesheet.

### Preparing samplesheet, gene reference and configuration file

To run the airrflow pipeline on single cell TCR or BCR sequencing data from fastq files, we need to prepare samplesheet, pre-built 10x genomics V(D)J references and configuration file in advance. Details on the required columns for this samplesheet are available [here](https://nf-co.re/airrflow/usage#fastq-input-samplesheet-single-cell-sequencing).

> [!WARNING]
> The fastq file names must follow the 10X Genomics file naming convention or the cellranger process will fail.

The prepared samplesheet for this tutorial is [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/10x_sc_raw.tsv) and a prepared configuration file is [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/resource.config). Download these two files to the directory where you intend to run the airrflow pipeline.

> [!TIP]
> Before setting memory and cpus in the configuration file, we recommend verifying the available memory and cpus on your system. Otherwise, exceeding the system's capacity may result in an error indicating that you requested more cpus than available or run out of memory.

Pre-built 10x genomics V(D)J references can be accessed at the [10x Genomics website](https://www.10xgenomics.com/support/software/cell-ranger/downloads). Both human and mouse V(D)J references are available. Download the reference that corresponds to the species of your dataset.

### Running airrflow

With all the files ready, it's time to run the airrflow pipeline.

```bash
nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode fastq \
--input 10x_sc_raw.tsv \
--library_generation_method sc_10x_genomics \
--reference_10x refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0 \
-c resource.config \
--clonal_threshold 0 \
--outdir sc_from_fastq_results \
-resume
```

In this tutorial, since the samples are TCRs, which do not have somatic hypermutation, clones are defined strictly by identical junction regions. For this reason, we set the `--clonal_threshold` parameter to 0. For more details on important considerations when performing clonal analysis check the section below.

Of course you can wrap all your code in a bash file. We prepared one for you and it's available [here](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/airrflow_sc_from_fastq.sh).
With the bash file, it's easy to run the pipeline with a single-line command.

```bash
bash airrflow_sc_from_fastq.sh
```

After launching the pipeline the following will be printed to the console output, followed by some Nextflow parameters and executions of Airrflow processes:

```bash
 N E X T F L O W   ~  version 24.10.5

Launching `https://github.com/nf-core/airrflow` [gloomy_monod] DSL2 - revision: d91dd840f4 [4.2.0]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/airrflow 4.2.0
------------------------------------------------------
```

Once the pipeline has finished successfully, the following message will appear:

```bash
-[nf-core/airrflow] Pipeline completed successfully-
Completed at: 11-Mar-2025 13:18:13
Duration    : 2m 46s
CPU hours   : 0.3 (0.1% cached)
Succeeded   : 17
Cached      : 2
```

## Important considerations for clonal analysis

An important step in the analysis of AIRR sequencing data is inferring B cell and T cell clones, or clonal groups, sometimes also called clonotypes. These are cells that are derived from the same progenitor cell through clonal expansion. For T cells, this definition is more strict as T cells do not undergo somatic hypermutation, so the TCRs from T cells in the same clone should be identical. For B cells, on the other hand, the BCRs from cells in the same clone can differ due to somatic hypermutation. They also can have a variety of isotypes.

There are two crucial considerations when defining clonal groups with nf-core/airrflow: across which samples should clonal groups be defined, and what should be the clonal threshold, i.e. how different can these receptors be, so that these are assigned to the same clonal group. These are discussed in detail in the following sections.

### Defining clonal groups across samples

Often times we want to analyze clonal groups from the same individual or animal model across time, different conditions or across samples extracted from different tissues. To ensure that the same clone ID (field `clone_id` in the output AIRR rearrangement file) is assigned to the same BCR / TCR clone across these conditions to be able to track the clones, the clonal inference step should be done pulling the sequences from these samples together. This is why, by default, nf-core/airrflow uses the `subject_id` column to group samples prior to defining clonal groups, so it is important to set the exact same subject ID to samples from the same individual across different conditions.

The sample grouping can also be controlled with the [`--cloneby`](https://nf-co.re/airrflow/4.2.0/parameters/#cloneby) parameter, by providing the name of the column containing the group information that should be used to pull the samples together before defining clonal groups (samples or rows with the same string in this column will be grouped together). You can create a new column if you wish for this purpose.

### Clonal inference method

nf-core/airrflow utilizes the Hierarchical clustering method in the [SCOPer](https://scoper.readthedocs.io/) Immcantation tool to infer clonal groups, which initially partitions the BCR / TCR sequences according to V gene, J gene and junction length. Then, it defines clonal groups within each partition by performing hierarchical clustering of the sequences within a partition and cutting the clusters according to an automatically detected or user-defined threshold. More details about this method can be found on the respective SCOPer [vignette](https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/#). Details on how to determine the clonal threshold can be found in the next section.

### Setting a clonal threshold

The clonal threshold can also be customized through the `--clonal_threshold` parameter. The clonal threshold specifies how different two BCRs can be so that are assigned to the same clonal group. The value is specified in length-normalized hamming distance across the BCR junction regions. By default, `--clonal_threshold` is set to be 'auto', allowing the clonal threshold to be determined automatically using a method included in the [SHazaM](https://shazam.readthedocs.io/) Immcantation tool. You can read more details about the method in the SHazaM [vignette](https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/).

For BCR data, we recommend using this default setting initially. After running the pipeline, you can review the automatically calculated threshold in the `find_threshold` report to make sure it is fitting the data appropriately. If the threshold is unsatisfactory, you can re-run the pipeline with a manually specified threshold (e.g. `--clonal_threshold 0.1`) that is appropriate for your data. For a low number of sequences that are insufficient to satisfactorily determine a threshold with this method, we generally recommend a threshold of 0.1 (length-normalized Hamming distance of nearest neighbors) for human BCR data.

Since TCRs do not undergo somatic hypermutation, TCR clones are defined strictly by identical junction regions. For this reason, the `--clonal_threshold` parameter should be set to 0 for TCR data.

### Including BCR lineage tree computation

BCR lineage tree computation is performed using the [Dowser](https://dowser.readthedocs.io/) Immcantation package. This step is skipped by default because it can be time-consuming depending on the size of the input data and the size of the clonal groups. To enable lineage tree computation, add the `--lineage_trees` parameter set to true. You can easily add lineage tree computation to a previous analysis by re-running the pipeline with the `-resume` so all the previous analysis steps are cached and not recomputed.

Dowser supports different methods for the lineage tree computation, `raxml` is the default but you can set other methods with the `--lineage_tree_builder` parameter, and provide the software executable with the `--lineage_tree_exec` parameter.

## Understanding the results

After running the pipeline, several subfolders are available under the results folder.

```bash
Airrflow_report.html
- cellranger
- vdj_annotation
- qc_filtering
- clonal_analysis
- repertoire_comparison
- multiqc
- report_file_size
- pipeline_info
```

The summary report, named `Airrflow_report.html`, provides an overview of the analysis results, such as an overview of the number of sequences per sample in each of the pipeline steps, the V(D)J gene assignment and QC, and V gene family usage. Additionally, it contains links to detailed reports for other specific analysis steps.

The analysis steps and their corresponding folders, where the results are stored, are briefly listed below. Detailed documentation on the pipeline output can be found on the [Output documentation page](https://nf-co.re/airrflow/docs/output/).

1. QC and sequence assembly (if starting from fastq files).

   - In this first step, Cell Ranger's VDJ algorithm is employed to assemble contigs, annotate contigs, call cells and generate clonoytpes. The results are stored in the 'cellranger' folder.

2. V(D)J annotation and filtering.

   - In this step, V(D)J gene segments are inferred using the provided germline reference and [`IgBLAST`](https://www.ncbi.nlm.nih.gov/igblast/). Alignments are annotated in AIRR format. Non-productive sequences and sequences with low alignment quality are filtered out unless otherwise specified. The intermediate results are stored under the folder named 'vdj_annotation'.

3. QC filtering.

   - In this step, cells without heavy chains or with multiple heavy chains are removed. Sequences in different samples that share the same cell_id and necleotide sequence are filtered out. The result are stored in the 'qc-filtering' folder.

4. Clonal analysis.

   - Results of the clonal threshold determination using `SHazaM` should be inspected in the html report under the 'clonal_analysis/find_threshold' folder. If the automatic threshold is unsatisfactory, you can set the threshold manually and re-run the pipeline.
     (Tip: use -resume whenever running the Nextflow pipeline to avoid duplicating previous work).
   - Clonal inference is performed with `SCOPer`. Clonal inference results as well as clonal abundance and diversity plots can be inspected in the html report in the folder 'clonal_analysis/define_clones'. For BCR sequencing data, mutation frequency is also computed using `SHazaM` at this step and plotted in the report. The `repertoires` subfolder contains the AIRR formatted files with the clonal assignments in a new column `clone_id` and mutation frequency in the column `mu_freq`. The `tables` subfolder contains the tabulated abundance and diversity computation as well as a table with the number of clones and their size. The `ggplots` subfolder contains the abundance and diversity plots as an `RData` object for loading and customization in R.
   - If lineage trees were computed using `Dowser`, a folder under 'clonal_analysis/dowser_lineages' will be present. The trees can be inspected in the html report and saved as PDF. Additionally, an `RDS` object with the formatted trees can also be loaded in R for customizing the lineage tree plots with Dowser.

5. Repertoire analysis.

   - Example calculation of several repertoire characteristics, e.g. V gene usage, for comparison between subjects, time points or cell populations is shown in the html report under `repertoire_comparison`. This report is generated from an Rmarkdown `Rmd` file. It is possible to customize this to fit the user's needs by editing the report and then providing the edited Rmd file with the `--report_rmd` parameter. Check also the remaining [Report parameters](https://nf-co.re/airrflow/parameters/#report-options) for further customizing the report.

6. Other reporting.
   Additional reports are also generated, including:
   - MultiQC report: summarizes QC metrics across all samples.
   - Pipeline_info report: various reports relevant to the running and execution of the pipeline.
   - Report_file_size report: Summary of the number of sequences left after each of the most important pipeline steps.

## Understanding error messages

Here, we list some common errors you may encounter while running the nf-core/airrflow pipeline and how to solve them.

1. Missing required column(s) in samplesheet.

   - The samplesheet collects experimental details that are important for the data analysis. Details on the required columns of a samplesheet are available [here](https://nf-co.re/airrflow/usage#assembled-input-samplesheet-bulk-or-single-cell-sequencing).

   - An example error message is shown below if the required column 'sex' is missing from the samplesheet( [assembled_samplesheet_missing_sex.tsv](https://github.com/nf-core/airrflow/blob/dev/docs/usage/single_cell_tutorial/sample_data_code/assembled_samplesheet_missing_sex.tsv) and the pipeline is run with this samplesheet.

```bash
#! /usr/bin/bash

nextflow run nf-core/airrflow -r 4.2.0 \
-profile docker \
--mode assembled \
--input assembled_samplesheet_missing_sex.tsv \
--outdir sc_from_assembled_results_error_test  \
-c resource.config \
-resume
```

```bash
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (assembled_samplesheet_missing_sex.tsv): Validation of file failed:
        -> Entry 1: Missing required field(s): sex
        -> Entry 2: Missing required field(s): sex
        -> Entry 3: Missing required field(s): sex
        -> Entry 4: Missing required field(s): sex

 -- Check script '/home/hl2244/.nextflow/assets/nf-core/airrflow/./workflows/../subworkflows/local/utils_nfcore_airrflow_pipeline/../../nf-core/utils_nfschema_plugin/main.nf' at line: 39 or see '.nextflow.log' file for more details
```

For more information on Nextflow errors and how to debug them you can check this [Nextflow troubleshooting tutorial](https://training.nextflow.io/2.1.1/basic_training/debugging/).

## Costumizing your analysis and generating your own figures

nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. The computations results (e.g. clonal inference, mutation frequency analysis) are stored in the output AIRR rearrangement repertoire files in newly generated columns under `clonal_analysis/define_clones/all_repertoires`. You can use these Airrflow results as input for customized analyses using R and the Immcantation tools. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html).

## Updating process resource requests

By default the pipeline has set reasonable process resource requests (number of CPUs, RAM memory, time limits) to the compute system. Depending on the size of your datasets or your running infrastructure you can customize these requests to meet your needs.

To update the resource requests for a specific pipeline process, you can do so in the `resource.config` file provided with the `-c` parameter. For example, to update the resource requests for the `CHANGEO_ASSIGNGENES` process:

```json title="resource.config"
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]

   withName:CHANGEO_ASSIGNGENES {
        cpus   = 2
        memory = 10.GB
        time   = 5h
   }
}
```

In nf-core pipelines, each process has a label indicating the resources that are being requested (`process_low`, `process_medium`, `process_high`, ...). The CPUs, RAM and time set up for each of these labels can be found in the [base.config](https://github.com/nf-core/airrflow/blob/master/conf/base.config) file. You can update the resource requests for all processes with a specific label by adding a new setting in your `resource.config` file provided with the `-c` parameter. For example here we update the process requests of processes with the `process_high` label:

```json title="resource.config"
process {
   resourceLimits = [cpus: 24, memory: 100.GB, time: 24.h]

   withLabel:process_high {
        cpus   = 24
        memory = 100.GB
        time   = 10h
   }
}
```

Note that the resource requests will never exceed what is specified in the `resourceLimits` line, so if you do want to increase the resource requests for specific processes, you should also increase the `resourceLimits` requests and run the pipeline in a compute infrastructure with sufficient resources. In this exmaple we also have updated the `resourceLimits` to reflect that.

> [!TIP]
> For more information about nf-core pipeline resource configurations, check out the [nf-core pipeline configuration docs](https://nf-co.re/docs/usage/getting_started/configuration).
