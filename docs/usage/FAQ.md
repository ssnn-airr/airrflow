# nf-core/airrflow: Frequently Asked Questions

## How to update process resource requests and resource limits?

By default, the pipeline defines reasonable resource requests for each process (number of CPUs, RAM memory, time limits) based on typical compute environments. However, you can adjust these settings to better match the size of your datasets or the capabilities of your compute infrastructure. You can customize the limits and requests in `resource.config` file and provide it to the pipeline using the -c parameter during execution. The `resourceLimits` option applies upper resource request limits to all the processes in the pipeline. Ensure that these limits do not exceed the available resources on your compute system.

```json title="resource.config"
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]
}
```

To update the resource requests for a specific pipeline process, you can also provide specific process requests in this config file. For example, to update the resource requests for the `CHANGEO_ASSIGNGENES` process:

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

In nf-core pipelines, each process has a label indicating the resources that are being requested (`process_low`, `process_medium`, `process_high`, ...). The CPUs, RAM and time set up for each of these labels can be found in the [base.config](https://github.com/nf-core/airrflow/blob/master/conf/base.config) file. You can update the resource requests for all processes with a specific label by providing the updated configuration. For example here we update the resource requests of processes with the `process_high` label:

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

Note that the resource requests will never exceed what is specified in the `resourceLimits` line, so if you do want to increase the resource requests for specific processes, you should also increase the `resourceLimits` requests and run the pipeline in a compute infrastructure with sufficient resources. In this example we also have updated the `resourceLimits` to reflect that.

> [!TIP]
> For more information about nf-core pipeline resource configurations, check out the [nf-core pipeline configuration docs](https://nf-co.re/docs/usage/getting_started/configuration).

## How to customize the analysis and figures?

nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. You can use nf-core/airrflow results as input for customized analyses using R and the Immcantation tools. There are three options to customize your analysis:

- Option 1: some of the intermediate analysis steps are stored on `RData` objects that can be loaded in R to customize your figures. For instance, clonal abundance calculations can be time-consuming, so the results are stored in the results folder (`clonal_abundance/define_clones/all_reps_clone_report/ggplots/abundanceSample.RData`). With `load()` function in R, both the abundance plot and the clonal abundance object can be loaded.
- Option 2: perform your own downstream analysis with the Immcantation framework. You can load the nf-core/airrflow results in AIRR format in R and use the Immcantation tools to plot the data as you need for publications. Check the [Immcantation tutorials](https://immcantation.readthedocs.io/en/stable/getting_started/getting-started.html) for this purpose, e.g. the Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html) shows an example for single-cell data analysis.
- Option 3: for more advanced users and in case you need to repeat the exact same analysis for multiple projects, you can customize the [Airrflow report](https://github.com/nf-core/airrflow/blob/master/assets/repertoire_comparison.Rmd) Rmarkdown file that comes with the pipeline and provide the updated version to the pipeline with the `--report_rmd` option. The pipeline will then use this file instead to create the report.
