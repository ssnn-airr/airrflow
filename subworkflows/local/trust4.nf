#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PARSE_TEMPLATE } from './modules/parse-template.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/fastqc.nf'
include { IMGT_TO_TRUST4 } from './modules/trust4.nf'
//include { IMGT_TRUST4_GENES } from './modules/trust4.nf'
//include { GENOME_TO_TRUST4 } from './modules/trust4.nf'
include { RUN_TRUST4 } from './modules/trust4.nf'
include { MULTI_TRUST4 } from './modules/trust4.nf'
include { CHANGEO_PARSEDB_SELECT } from './modules/trust4.nf'
include { AIRR_SAMPLESHEET as IG_AIRR_SAMPLESHEET } from './modules/trust4.nf'
include { AIRR_SAMPLESHEET as TR_AIRR_SAMPLESHEET } from './modules/trust4.nf'
include { MERGE_SAMPLESHEETS as IG_MERGE_SAMPLESHEETS } from './modules/trust4.nf'
include { MERGE_SAMPLESHEETS as TR_MERGE_SAMPLESHEETS } from './modules/trust4.nf'

log.info """\
Run cellranger vdj
===================================
Project parameters:
- samplesheet                : ${params.samplesheet}
- organism                   : ${params.organism}
- outdir                     : ${params.outdir}
"""

workflow {

    // Create channel from input samplesheet
    Channel.fromPath("${params.samplesheet}")
            .ifEmpty { exit 1, "File not foud: ${params.samplesheet}" }
            .set { ch_samplesheet }
    //ch_samplesheet
    //.dump(tag: 'ch_samplesheet')

    /* Not used
    // Create channels for reference genomes
    //if( params.genome ){
    //  Channel.fromPath("${params.genome}")
    //         .ifEmpty { exit 1, "Genome not found: ${params.genome}" }
    //         .set { ch_genome }
    //} else {
    //  ch_genome = Channel.value()
    //}
    //ch_genome
    //.dump(tag: 'ch_genome')

    //if( params.genome ){
    //  Channel.fromPath("${params.genome_annotation}")
    //         .ifEmpty { exit 1, "Genome annotation not found: ${params.genome_annotation}" }
    //         .set { ch_genome_annotation }
    //} else {
    //  ch_genome_annotation = Channel.value()
    //}
    //ch_genome_annotation
    //.dump(tag: 'ch_genome_annotation')
    */

    PARSE_TEMPLATE(ch_samplesheet)
    //PARSE_TEMPLATE.out.ch_template
    //.dump(tag: 'parsed template')


    ch_samples = PARSE_TEMPLATE.out.ch_template
    .splitCsv(header:true, sep:'\t', quote:'"')
    .map{ row-> [
                    meta:[
                        'sample_id':row.'sample_name',
			'subject_id':row.'subject_name',
			'organism':row.'organism',
			'tissue':row.'tissue',
			'biological_sex':row.'biological_sex',
			'age':row.'age',
			'collaborator':row.'collaborator',
                        'sample_id_alt': row.'sample_alternate_name'.tokenize( ',' )[0],
                        'experiment_type':row.'experiment_type'],
                    reads:[
                        'R1':row.R1 != 'NA' ? file(row.R1.trim()) : "",
                        'R2':row.R2 != 'NA' ? file(row.R2.trim()) : "",
                        'R3':""
                        ]
                ]}

    //ch_samples
    //.dump(tag: 'ch_samples')

    // Run fastQC
    FASTQC(ch_samples)
    MULTIQC(FASTQC.out.zip.collect())

    // Create TRUST4 reference germline and genome database
    IMGT_TO_TRUST4("${params.organism}")

    // Not used
    //IMGT_TRUST4_GENES(IMGT_TO_TRUST4.out.imgt_trust4)
    //GENOME_TO_TRUST4(ch_genome.collect(),
    //                  ch_genome_annotation.collect(),
    //                  IMGT_TRUST4_GENES.out.trust4_genes)
    //

    RUN_TRUST4(
      ch_samples,
      IMGT_TO_TRUST4.out.imgt_trust4.collect() // ,GENOME_TO_TRUST4.out.trust4_genome.collect()
    )

    MULTI_TRUST4(
        RUN_TRUST4.out.tab
        .map{ it -> [ it[1] ] } // *_airr.tsv
        .collect(),
        RUN_TRUST4.out.trust4_results
        .map{ it -> [ it[1] ] } // *_read_count.txt
        .collect(),
        RUN_TRUST4.out.tab
        .map{ it -> [ it[0]['sample_id'] ] } // sample_id
        .collect()
    )

    CHANGEO_PARSEDB_SELECT(
        RUN_TRUST4.out.tab
    )

    CHANGEO_PARSEDB_SELECT.out.ig_tab
    .dump(tag:'ig_tab')

    IG_AIRR_SAMPLESHEET(
        CHANGEO_PARSEDB_SELECT.out.ig_tab,
        "ig"  // receptor: ig or tr
    )
    TR_AIRR_SAMPLESHEET(
        CHANGEO_PARSEDB_SELECT.out.tr_tab,
        "tr" // receptor: ig or tr
    )

    IG_MERGE_SAMPLESHEETS(IG_AIRR_SAMPLESHEET.out.samplesheet.collect{it[1]}, "ig")
    TR_MERGE_SAMPLESHEETS(TR_AIRR_SAMPLESHEET.out.samplesheet.collect{it[1]}, "tr")
}
