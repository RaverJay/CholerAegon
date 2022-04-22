nextflow.enable.dsl=2


//// TODO
// - example testdata run to make a 'local install'


// defaults, these are overridden by the 'nextflow run' invocation on the command line
params.fasta = false
params.samples = false
params.longreads = false
params.shortreads = false

params.do_all_assemblies = false

params.min_coverage_percent = 80
params.min_identity_percent = 80


// channels

params.rgi_card_localDB = workflow.projectDir + "/data/CARD/localDB"
rgi_card_localDB_ch = Channel.fromPath(params.rgi_card_localDB , checkIfExists: true)

params.card_ontology_obo = workflow.projectDir + "/data/CARD/aro.obo"
card_ontology_obo_ch = Channel.fromPath(params.card_ontology_obo , checkIfExists: true)


if (params.genome_reference) {
    genome_reference_ch = Channel
        .fromPath(params.genome_reference , checkIfExists: true)
}


///////////////////
//// Processes ////
///////////////////

////
// preprocess

process concat_fastq {

    label 'base'
    cpus = 1
    input:
        tuple val(name), val(fastqs)
    output:
        tuple val(name), path("${name}_concat.fastq")
    shell:
        '''
        NUM=`wc -w <<< "!{fastqs}"`
        echo $NUM
        if [[ $NUM -ge 2 ]]; then
        cat !{fastqs} > !{name}_concat.fastq
        else
        ln -s !{fastqs} !{name}_concat.fastq
        fi
        '''
}


process prepare_fasta {

    label 'base'
    cpus = 1
    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path("fixed_names/${name}.fasta")
    shell:
        '''
        mkdir fixed_name
        for F in *.fa *.fasta *.fna; do
        if [ -e "$F" ]; then
        ln -s $F fixed_names/!{name}.fasta
        '''
}



////
// longreads

process flye {

    label 'flye'
    cpus = "${params.cpus}"
    memory = "${params.memory}"
    publishDir "${params.output}/sample_results/${name}/${type}", mode: 'copy', pattern: 'flye_assembly'
    errorStrategy 'ignore'

    input:
        tuple val(name), val(type), path(reads)
    output:
        tuple val(name), val(type), path('flye_assembly')
        tuple val(name), val(type), path("${name}_flye.fasta"), emit: assembly

    script:
        """
        flye --nano-raw ${reads} --plasmids --out-dir flye_assembly --threads ${params.cpus}
        cp flye_assembly/assembly.fasta "${name}_flye.fasta"
        """
}


process minimap2_paf {

    label 'minimap2'
    cpus = "${params.cpus}"

    input:
        tuple val(name), path(reads), path(target)
    output:
        tuple val(name), path("${name}.paf")

    script:
        """
        minimap2 -t ${params.cpus} -x map-ont ${target} ${reads} > ${name}.paf
        """
}


process medaka {

    label 'medaka'
    cpus = "${params.cpus}"
    publishDir "${params.output}/sample_results/${name}/${type}", mode: 'copy', pattern: 'medaka_consensus'
    publishDir "${params.output}/sample_results/${name}/", mode: 'copy', pattern: "${name}_${type}_assembly.fasta"

    input:
        tuple val(name), val(type), path(reads), path(draft_assembly)
    output:
        tuple val(name), val(type), path('medaka_consensus')
        tuple val(name), val(type), path("${name}_${type}_assembly.fasta"), emit: assembly
    script:
        // medaka will never use more than 2 cpu cpus - TODO use GPU if available
        """
        medaka_consensus -i ${reads} -d ${draft_assembly} -o medaka_consensus -t ${params.cpus} -m ${params.medaka_model}
        cp medaka_consensus/consensus.fasta ${name}_${type}_assembly.fasta

        # cleanup to save disk space
        rm medaka_consensus/consensus_probs.hdf
        rm medaka_consensus/calls_to_draft.bam
        rm medaka_consensus/calls_to_draft.bam.bai
        """
}


////
// shortreads


process fastp {

    label 'fastp'
    cpus = "${params.cpus}"
    memory = "${params.memory}"
    publishDir "${params.output}/sample_results/${name}", mode: 'copy', pattern: 'fastp_log'

    input:
        tuple val(name), path(shortreads_p1), path(shortreads_p2)
    output:
        tuple val(name), path("${name}_p1_fastp.fq.gz"), path("${name}_p2_fastp.fq.gz")

    script:
        """
        mkdir fastp_log
        fastp -5 -3 -W 3 -M 20 -l 30 --in1 ${shortreads_p1} --in2 ${shortreads_p2} \
            --out1 ${name}_p1_fastp.fq.gz --out2 ${name}_p2_fastp.fq.gz \
            -h fastp_log/${name}_fastp_log.html
        """
}


process spades_denovo {

    label 'spades'
    cpus = "${params.cpus}"
    memory = "${params.memory}"
    publishDir "${params.output}/sample_results/${name}/${type}", mode: 'copy', pattern: 'spades_assembly'
    publishDir "${params.output}/sample_results/${name}/", mode: 'copy', pattern: "${name}_${type}_assembly.fasta"

    input:
        tuple val(name), val(type), path(shortreads_p1), path(shortreads_p2)
    output:
        tuple val(name), val(type), path("spades_assembly")
        tuple val(name), val(type), path("${name}_${type}_assembly.fasta"), emit: assembly
    script:
        """
        spades.py -t ${task.cpus} --isolate -1 ${shortreads_p1} -2 ${shortreads_p2} -o spades_assembly
        cp spades_assembly/scaffolds.fasta ${name}_${type}_assembly.fasta
        """
}


////
// hybrid


process bwamem2_paired {

    label 'bwamem2'
    cpus = "${params.cpus}"

    input:
        tuple val(name), path(illumina_p1), path(illumina_p2), path(assembly)
    output:
        tuple val(name), path("${name}_bwamem_vs_ref.sam")
    script:
        """
        bwa-mem2 index ${assembly}
        bwa-mem2 mem -t ${params.cpus} -o ${name}_bwamem_vs_ref.sam ${assembly} ${illumina_p1} ${illumina_p2} 
        """
}


process samtools_bam {

    label 'samtools'
    cpus = 1

    input:
        tuple val(name), path(samfile)
    output:
        tuple val(name), path("*.bam"), path("*.bam.bai")
    shell:
        '''
        set -eo pipefail

        SAMIN=!{samfile}
        samtools view -hbF4 $SAMIN | samtools sort > ${SAMIN%.sam}.bam && rm $SAMIN
        samtools index ${SAMIN%.sam}.bam
        '''
}


process pilon {

    label 'pilon'
    cpus = "${params.cpus}"
    memory = '8 GB'
    publishDir "${params.output}/sample_results/${name}/${type}", mode: 'copy', pattern: 'pilon_polishing'
    publishDir "${params.output}/sample_results/${name}/", mode: 'copy', pattern: "${name}_${type}_assembly.fasta"

    input:
        tuple val(name), val(type), path(mapping), path(mapping_idx), path(assembly)
    output:
        tuple val(name), val(type), path("pilon_polishing")
        tuple val(name), val(type), path("${name}_${type}_assembly.fasta"), emit: assembly
    script:
        // NEEDS 8G memory machine for now
        """
        pilon -Xmx8g --threads "${params.cpus}" --genome ${assembly} --frags ${mapping} --outdir pilon_polishing
        cp pilon_polishing/pilon.fasta ${name}_${type}_assembly.fasta
        """
}


////
// AMR gene detection


process resistance_gene_identifier {

    label 'rgi'
    cpus = "${params.cpus}"
    publishDir "${params.output}/sample_results/${name}/AMR/", mode: 'copy'

    input:
        tuple val(name), val(type), path(assembly), path(localDB)
    output:
        tuple val(name), val(type), path("RGI_${name}_${type}.txt")
    script:
        """
        rgi main -n ${params.cpus} --input_sequence ${assembly} --output_file RGI_${name}_${type} --input_type contig --local
        """
}

process aggregate_rgi_results {

    label 'base'
    cpus = 1
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(rgi_results)
    output:
        path("RGI_aggregated_results.csv")
    script:
        """
        aggregate_rgi_results.py ${rgi_results} > RGI_aggregated_results.csv
        """
}


process abricate {

    label 'abricate'
    cpus = "${params.cpus}"
    publishDir "${params.output}/sample_results/${name}/AMR/", mode: 'copy'

    input:
        tuple val(name), val(type), path(assembly)
    output:
        tuple val(name), val(type), path("Abricate_CARD_${name}_${type}.csv")
    script:
        """
        abricate --threads ${params.cpus} --db card ${assembly} > Abricate_CARD_${name}_${type}.csv
        """
}


process aggregate_abricate_results {

    label 'abricate'
    cpus = "${params.cpus}"
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(abricate_results)
    output:
        path("Abricate_CARD_aggregated_results.csv")
    script:
        """
        abricate --threads ${params.cpus} --summary ${abricate_results} > Abricate_CARD_aggregated_results.csv
        """
}


process combine_results {

    label 'base'
    cpus = 1
    publishDir "${params.output}/sample_results/${name}/AMR/", mode: 'copy'

    input:
        tuple val(name), val(type), path(rgi_results), path(abricate_results)
    output:
        tuple val(name), val(type), path("${name}_${type}_combined_results.csv")
    script:
        """
        combine_results.py --rgi_results ${rgi_results} \
            --abricate_results ${abricate_results} \
            --output ${name}_${type}_combined_results \
            --coverage_threshold ${params.min_coverage_percent} \
            --identity_threshold ${params.min_identity_percent}
        """
}


process aggregate_combined_results {

    label 'base'
    cpus = 1
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(combined_results)
    output:
        path("AMR_combined_aggregated_results.csv")
    script:
        """
        aggregate_combined_results.py ${combined_results} > AMR_combined_aggregated_results.csv
        """
}


process predict_resistances {

    label 'base'
    cpus = 1
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(combined_agg_amr_results)
        path(card_ontology_obo)
    output:
        path("drug_resistance_prediction.tsv")
    script:
        """
        predict_drug_resistance.py --aro_ontology ${card_ontology_obo} --gene_results ${combined_agg_amr_results} --output drug_resistance_prediction.tsv
        """
}



////
// stats and summary


process fastani {

    label 'fastani'
    cpus = "${params.cpus}"
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(assemblies)
        path(genome_reference)
    output:
        path("fastANI_aggregated_results.tsv")
    script:
        """
        ls *.fasta > assemblies_list.txt
        fastANI --threads ${params.cpus} -r ${genome_reference} --ql assemblies_list.txt -o fastANI_aggregated_results.tsv
        """
}

process summary_report {

    label 'base'
    cpus = 2
    publishDir "${params.output}/", mode: 'copy'

    input:
        path(amr_results)
        path(fastani_results)
    output:
        path("CholerAegon_summary_report_*.html")
    script:
        """
        summary_report.py --amr_results ${amr_results} --fastani_results ${fastani_results}
        """
}


///////////////
// Workflows //
///////////////


workflow longread_assembly_wf {

    take:
        longreads

    main:
        // flye denovo
        input_lreads = longreads.map{ it -> [ it[0], 'longreads', it[1] ] }
        flye( input_lreads )

        // medaka consensus polish
        medaka(input_lreads.join( flye.out.assembly.map{ it -> [ it[0], it[2] ] } ) )

    emit:
        medaka.out.assembly
}


workflow shortread_assembly_wf {

    take:
        shortreads

    main:
        // spades assembly
        input_sreads = shortreads.map{ it -> [ it[0], 'shortreads', it[1], it[2] ] }
        spades_denovo( input_sreads )

    emit:
        spades_denovo.out.assembly

}


workflow hybrid_assembly_wf {

    take:
        longreads_assembly
        shortreads

    main:
        // pilon hybrid polish
        input_bwa = shortreads.join( longreads_assembly.map{ it -> [ it[0], it[2] ] } )
        input_pilon = samtools_bam( bwamem2_paired( input_bwa ) ).join( longreads_assembly.map{ it -> [ it[0], it[2] ] } )
        input_pilon_type = input_pilon.map{ it -> [ it[0], 'hybrid', it[1], it[2], it[3] ] }
        pilon( input_pilon_type )

    emit:
        pilon.out.assembly
}


workflow res_gene_detection_wf {

    take:
        assemblies

    main:
        // RGI
        resistance_gene_identifier( assemblies.combine( rgi_card_localDB_ch ) )
        aggregate_rgi_results( resistance_gene_identifier.out.map{ it -> it[2] }.collect() )
        
        // Abricate
        abricate( assemblies )
        aggregate_abricate_results( abricate.out.map{ it -> it[2] }.collect() )

        // Combined
        combine_results( resistance_gene_identifier.out.join( abricate.out.map{ it -> [ it[0], it[2] ] } ) )
        aggregate_combined_results( combine_results.out.map{ it -> it[2] }.collect() )

    emit:
        aggregate_combined_results.out

}


workflow res_prediction_wf {

    take:
        combined_agg_amr_results
        card_ontology

    main:
        predict_resistances( combined_agg_amr_results, card_ontology )

    emit:
        predict_resistances.out

}


workflow create_summary_report_wf {

    take:
        predicted_resistances
        fastani_results

    main:
        summary_report( predicted_resistances, fastani_results )

}


//////////
// Main //
//////////

workflow {
    main:

        //////////////
        // Assembly //
        //////////////

        // fasta input
        if (params.fasta) {
            fasta_input_ch = Channel
                .fromPath( params.fasta, checkIfExists: true)
                .map{ file -> tuple( file.simpleName, file ) }

            assemblies_ch = prepare_fasta( fasta_input_ch )
        }

        else {

            // samples input
            if (params.samples) {
                samples_input_ch = Channel
                    .fromPath( params.samples, checkIfExists: true )
                    .splitCsv()
                    .map { row -> ["${row[0]}", "${row[1]}",
                                file("${row[2]}", checkIfExists: true), file("${row[3]}", checkIfExists: true)] }
                    // csv table with columns:  samplename, longreads, illumina_p1, illumina_p2

                longread_input_ch = concat_fastq( samples_input_ch.map{ it -> [it[0], it[1]] } )
                shortread_input_ch = fastp( samples_input_ch.map{ it -> [it[0], it[2], it[3]] } )

                longread_assembly_wf( longread_input_ch )
                hybrid_assembly_wf( longread_assembly_wf.out, shortread_input_ch )

                if (! params.do_all_assemblies) {
                    assemblies_ch = hybrid_assembly_wf.out
                }
                else {
                    shortread_assembly_wf( shortread_input_ch )
                    assemblies_ch = hybrid_assembly_wf.out.mix( longread_assembly_wf.out ).mix( shortread_assembly_wf.out )
                }

                
            }

            else {
                // longreads input
                if (params.longreads) {
                    longread_input_ch = Channel
                        .fromPath( params.longreads, checkIfExists: true)
                        .map { file -> tuple(file.simpleName, file) }
                        // samplename, longreads

                    longread_assembly_wf( longread_input_ch )
                }

                // shortreads input
                // TODO unpaired input
                if (params.shortreads) {
                    shortread_input_raw_ch = Channel
                        .fromPath( params.longreads, checkIfExists: true)
                        .map { file -> tuple(file.simpleName, file) }
                        // samplename, shortreads

                    shortread_input_ch = fastp ( shortread_input_raw_ch )
                }

            
                if (params.longreads &&! params.shortreads) {
                    assemblies_ch = longread_assembly_wf.out
                }
                if (!params.longreads && params.shortreads) {
                    shortread_assembly_wf( shortread_input_ch )
                    assemblies_ch = shortread_assembly_wf.out
                }
                if (params.longreads && params.shortreads) {
                    hybrid_assembly_wf( longread_input_ch, shortread_input_ch )
                    assemblies_ch = hybrid_assembly_wf.out
                }
            }
        }
        // TODO add Guppy container + fast5 input

        
        ///////////////////
        // AMR detection //
        ///////////////////

        // detect genes
        res_gene_detection_wf( assemblies_ch )

        // predict drug resistances
        res_prediction_wf( res_gene_detection_wf.out, card_ontology_obo_ch )


        /////////////
        // Summary //
        /////////////

        if (params.genome_reference) {
            ani_result_ch = fastani( assemblies_ch.map{ it -> it[2] }.collect(), genome_reference_ch )
        } else {
            ani_result_ch = Channel.from( ['deactivated'] )
        }

        create_summary_report_wf( res_prediction_wf.out, ani_result_ch )


}
