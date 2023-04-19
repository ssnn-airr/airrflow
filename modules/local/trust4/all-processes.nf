process IMGT_TO_TRUST4 {
    tag 'imgt_to_trust4'
    label 'immcantation'

    publishDir "$params.outdir/imgt_trust4_ref", mode:"copy"

    input:
    val organism

    output:
    path("*.fa"), emit: imgt_trust4
    path("germlines"), emit: ref_germlines
    path("igblast"), emit: ref_igblast

    script:
    """
    #!/usr/bin/env Rscript

    imgtConvert <- function(..., from="immcantation", to="trust4") {
    inputs <- c(...)
    if (from=="immcantation") {
          if (to=="trust4") {
             Reduce(c, sapply(inputs, tigger::readIgFasta, USE.NAMES = T))
          } else {
             stop("Unknown to")
          }
       } else {
          stop("Unknown from")
       }
    }
    remotes::install_github("HenrikBengtsson/R.utils", ref="master", upgrade="never")
    library(R.utils)
    copyDirectory("/usr/local/share/germlines", "./germlines")
    copyDirectory("/usr/local/share/igblast", "./igblast")

    imgt_files <- c(
      list.files(file.path("/usr/local/share/germlines/imgt","$organism","constant"), full.names = T),
      list.files(file.path("/usr/local/share/germlines/imgt","$organism","vdj"), full.names = T)
    )
    trust4_imgt <- imgtConvert(imgt_files,
                               from="immcantation",
                               to="trust4")
    trust4_imgt_fn <- file.path(getwd(),"immcantation_IMGT+C.fa")
    tigger::writeFasta(trust4_imgt, file=trust4_imgt_fn)
    """
}

process IMGT_TRUST4_GENES {
    tag 'imgt_trust4_genes'
    label 'immcantation'

    publishDir "$params.outdir/imgt_trust4_ref", mode:"copy"

    input:
    path imgt_trust4

    output:
    path("*gene_name.txt"), emit: trust4_genes

    script:
    """
    grep ">" "${imgt_trust4}" | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt
    """
}

process GENOME_TO_TRUST4 {
    tag 'genome_to_trust4'
    label 'trust4'

    publishDir "$params.outdir/imgt_trust4_ref", mode:"copy"

    input:
    path genome
    path genome_annotation
    path genes

    output:
    path("bcrtcr.fa"), emit: trust4_genome
    tuple path("trust4_genome.fa"), path("trust4_genes.gtf")

    script:
    """
    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < "${genome}" \
    | grep "^>chr[1-9XY]*\$" -A 1 > trust4_genome.fa && \
    grep "^chr[1-9XY]\\+\\s" "${genome_annotation}" > trust4_genes.gtf && \
    perl /usr/local/bin/BuildDatabaseFa.pl trust4_genome.fa trust4_genes.gtf "${genes}" > bcrtcr.fa
    """
}

process RUN_TRUST4 {
    tag "${meta.'sample_id'}"
    label 'trust4'
    label 'process_high'

    publishDir "$params.outdir/trust4_out", mode:"copy"

    input:
    tuple val(meta), val(reads)
    path imgt_ref

    output:
    tuple val(meta), path("*_airr.tsv"), emit: tab
    tuple path("*_report.tsv"), path("*_read-count.txt"), path("*_raw.out"), path("*_final.out"), path("*_annot.fa"),path("*_cdr3.out"), path("*_assembled_reads.fa"), emit: trust4_results

    script:
    """
    run-trust4 --ref "${imgt_ref}" -f "${imgt_ref}" \
    -1 "${reads.R1}" \
    -2 "${reads.R2}" \
    -o "${meta.'sample_id'}" \
    -t ${task.cpus}
    if [[ "${reads.R1}" =~ \\.gz\$ ]]
    then
    zcat "${reads.R1}" | echo "R1 reads " \$((`wc -l`/4)) > "${meta.'sample_id'}_read-count.txt"
    else
    cat "${reads.R1}" | echo "R1 reads " \$((`wc -l`/4)) > "${meta.'sample_id'}_read-count.txt"
    fi
    if [[ "${reads.R2}" =~ \\.gz\$ ]]
    then
    zcat "${reads.R2}" | echo "R2 reads " \$((`wc -l`/4)) >> "${meta.'sample_id'}_read-count.txt"
    else
    cat "${reads.R2}" | echo "R2 reads " \$((`wc -l`/4)) >> "${meta.'sample_id'}_read-count.txt"
    fi
    grep "Found [0-9]* reads"  .command.log >> "${meta.'sample_id'}_read-count.txt"
    grep "Assembled [0-9]* reads"  .command.log >> "${meta.'sample_id'}_read-count.txt"
    grep "Try to rescue [0-9]* reads"  .command.log >> "${meta.'sample_id'}_read-count.txt"
    grep "Rescued [0-9]* reads"  .command.log >> "${meta.'sample_id'}_read-count.txt"
    """
}


process MULTI_TRUST4 {
    tag 'multi_trust4'
    label 'immcantation'

    publishDir "$params.outdir/trust4_out/multi_trust4", mode:"copy"

    input:
    path airr
    path read_count
    val sample_id

    output:
    path("*.html")

    script:
    """
    echo "${airr.join('\n')}" > input_airr.txt
    echo "${read_count.join('\n')}" > input_read_count.txt
    echo "${sample_id.join('\n')}" > input_sample_id.txt
    Rscript -e "wd<-getwd(); \\
        rmarkdown::render(\\
            '${projectDir}/assets/multi_trust4.Rmd', \\
            output_file = 'multi_trust4.html', \\
            knit_root_dir = wd, \\
            output_dir = wd, \\
            params = list( \\
                airr = 'input_airr.txt', \\
                read_count='input_read_count.txt', \\
                sample_id = 'input_sample_id.txt' ) \\
        )"
    """
}

process CHANGEO_PARSEDB_SELECT {
    tag "$meta.sample_id"
    label 'immcantation'
    publishDir "$params.outdir/trust4_out/parsedb-select", mode:"copy"

    input:
    tuple val(meta), path(tab) // sequence tsv in AIRR format

    output:
    tuple val(meta), path("*IG_parse-select.tsv"), emit: ig_tab, optional: true // sequence tsv in AIRR format
    tuple val(meta), path("*TR_parse-select.tsv"), emit: tr_tab, optional: true // sequence tsv in AIRR format
    path("*_command_log.txt"), emit: logs, optional: true //process logs

    script:
    """
    { # try
        ParseDb.py select -d "${tab}" -f v_call -u IG --regex --outname "${meta.sample_id}_IG" > "${meta.sample_id}_IG_select_command_log.txt"
    } || { # catch
        echo "Emplty file ${tab}" > "${meta.sample_id}_IG_select_command_log.txt" && rm "${meta.sample_id}_IG_parse-select.tsv"
    }
    { # try
        ParseDb.py select -d "${tab}" -f v_call -u TR --regex --outname "${meta.sample_id}_TR" > "${meta.sample_id}_TR_select_command_log.txt"
    } || { # catch
        echo "Emplty file ${tab}" > "${meta.sample_id}_TR_select_command_log.txt" && rm "${meta.sample_id}_TR_parse-select.tsv"
    }
    """
}


process AIRR_SAMPLESHEET {

    tag "${meta.'sample_id'}"
    label 'immcantation'
    //publishDir "$params.outdir/trust4_out/parsedb-select/", mode:"copy"

    input:
    tuple val(meta), val(airr)
    val(receptor)

    output:
    tuple val(meta), path("*_samplesheet.tsv"), emit: samplesheet

    exec:
    def fields = [
        filename : "$params.outdir/trust4_out/parsedb-select/${airr.getName()}",
        species : "${meta.organism}",
	    subject_id : "${meta.'subject_id'}",
	    sample_id : "${meta.'sample_id'}",
	    tissue : "${meta.'tissue'}",
	    sex : "${meta.'biological_sex'}",
	    age : "${meta.'age'}",
	    biomaterial_provider : "${meta.'collaborator'}",
	    pcr_target_locus : "${receptor}",
	    single_cell : "false"
    ]

    samplesheet  = fields.keySet().collect{it}.join("\t") + '\n'
    samplesheet += fields.values().collect{it}.join("\t")

    def samplesheet_file = task.workDir.resolve("${meta.'sample_id'}_${receptor}_samplesheet.tsv")
    samplesheet_file.text = samplesheet
}

process MERGE_SAMPLESHEETS {
    tag "all-samples"
    label 'immcantation'
    publishDir "$params.outdir/trust4_out/", mode:"copy"

    input:
    path ('samplesheets/*')
    val(receptor)

    output:
    path "*airrflow-samplesheet.tsv", emit: samplesheet

    script:
    """
    head -n 1 `ls ./samplesheets/* | head -n 1` > "${receptor}_airrflow-samplesheet.tsv"
    for fileid in `ls ./samplesheets/*`; do
        awk 'NR>1' \$fileid >> "${receptor}_airrflow-samplesheet.tsv"
    done
    """
}
