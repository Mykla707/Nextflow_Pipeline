#!/usr/bin/env nextflow

params.data_dir            = "/mnt/c/Users/mylen/Documents/BIOINF_BIO2025/Nextflow/nf-scripts/nextflow_pipeline/exam/exam/hepatitis"
params.glob                = "*.{fa,fasta}"
params.output_combined_file = "comb.fasta"
params.accession_id        = "M21012"
params.output_mafft        = "mafft.fasta"
params.output_trimal       = "trimal_.fasta"
params.trimal_rep_html     = "trimal_rep.html"
params.outdir              = "results"



// Process: Extract FASTA from NCBI
process fasta_data {
    conda 'bioconda::entrez-direct=24.0'

    input:
    val accession_id

    output:
    path "${accession_id}.fasta", emit: sequence_fasta

    script:
    """
    esearch -db nucleotide -query "${accession_id}" \
        | efetch -format fasta \
        > "${accession_id}.fasta"
    """
}

// Process: Combine input FASTA files
process comb_seq {
    
    input:
    path fasta_files

    output:
    path "${params.output_combined_file}", emit: combined_output

    script:
    """
    cat ${fasta_files.join(' ')} > ${params.output_combined_file}
    """
}

// Process: MAFFT alignment
process mafft_align {
    conda 'bioconda::mafft=7.525'

    input:
    path input_fasta

    output:
    path "${params.output_mafft}", emit: align_fasta

    script:
    """
    mafft --auto --thread -1 "${input_fasta}" > "${params.output_mafft}"
    """
}

// Process: Trimal trimming

process trimal_align {
    conda 'bioconda::trimal=1.5.0'

    publishDir "${params.outdir}/trimal_results", mode: 'copy', pattern: '*'

    input:
    path input_align_fasta

    output:
    path "${params.output_trimal}", emit: trimm_fasta
    path "${params.trimal_rep_html}", emit: trimal_report

    script:
    """
    trimal -in "${input_align_fasta}" -out "${params.output_trimal}" -automated1 -htmlout "${params.trimal_rep_html}"
    """
}

// MAIN WORKFLOW
workflow {

    // Download a FASTA file from NCBI
    fasta_result = fasta_data(params.accession_id)

    // Get all .fasta or .fa files from input directory
    input_files_ch = Channel.fromPath("${params.data_dir}/${params.glob}")

    // Combine downloaded file + local files into one channel
    all_fasta_ch = input_files_ch.mix(fasta_result.sequence_fasta)

    // Collect into a list of files for comb_seq
    combined_input_ch = all_fasta_ch.collect()

    // Combine sequences into one file
    combined = comb_seq(combined_input_ch)

    // Align combined fasta
    aligned = mafft_align(combined.combined_output)

    // Trim alignment and generate report
    trimmed = trimal_align(aligned.align_fasta)

}
