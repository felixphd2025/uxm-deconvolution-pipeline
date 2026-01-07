#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
    ========================================
    UXM DECONVOLUTION PIPELINE
    ========================================
    Samples         : ${params.sample_sheet}
    Atlas U25       : ${params.atlas_u25}
    Atlas U250      : ${params.atlas_u250}
    Output directory: ${params.outdir}
    ========================================
    """
    .stripIndent()

/*
 * PROCESS 1: Merge BAMs (if needed)
 */
process MERGE_BAMS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/merged", mode: 'copy'
    
    cpus 32
    memory '128 GB'
    
    input:
    tuple val(sample_id), path(bams)
    
    output:
    tuple val(sample_id), path("${sample_id}_merged.bam"), path("${sample_id}_merged.bam.bai"), emit: merged
    path "${sample_id}_merge.log", emit: log
    
    script:
    """
    echo "========================================" | tee ${sample_id}_merge.log
    echo "Merging ${sample_id}" | tee -a ${sample_id}_merge.log
    echo "Total BAM files: \$(ls -1 *.bam | wc -l)" | tee -a ${sample_id}_merge.log
    echo "Started: \$(date)" | tee -a ${sample_id}_merge.log
    echo "========================================" | tee -a ${sample_id}_merge.log
    
    # Create BAM list
    ls -1 *.bam > bam_list.txt
    
    # Merge
    samtools merge -@ ${task.cpus} -b bam_list.txt ${sample_id}_merged.bam
    
    # Index
    samtools index -@ ${task.cpus} ${sample_id}_merged.bam
    
    echo "Completed: \$(date)" | tee -a ${sample_id}_merge.log
    echo "Size: \$(ls -lh ${sample_id}_merged.bam | awk '{print \$5}')" | tee -a ${sample_id}_merge.log
    """
}

/*
 * PROCESS 2: BAM to PAT
 */
process BAM_TO_PAT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/pat", mode: 'copy'
    
    cpus 16
    memory '64 GB'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("*.pat.gz"), emit: pat
    path "${sample_id}_bam2pat.log", emit: log
    
    script:
    """
    export PATH=\$HOME/bioinformatics/wgbs_tools:\$PATH
    
    echo "Converting BAM to PAT for ${sample_id}" | tee ${sample_id}_bam2pat.log
    echo "Started: \$(date)" | tee -a ${sample_id}_bam2pat.log
    
    wgbstools bam2pat \
        ${bam} \
        --genome hg38 \
        --out_dir . \
        --nanopore \
        --np_thresh 0.8 \
        2>&1 | tee -a ${sample_id}_bam2pat.log
    
    echo "Completed: \$(date)" | tee -a ${sample_id}_bam2pat.log
    """
}

/*
 * PROCESS 3: UXM Deconvolution - U25 Atlas
 */
process UXM_DECONV_U25 {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/deconvolution", mode: 'copy'
    
    cpus 4
    memory '16 GB'
    
    input:
    tuple val(sample_id), path(pat)
    path atlas
    
    output:
    tuple val(sample_id), path("${sample_id}_deconv_U25.csv"), emit: deconv
    path "${sample_id}_deconv_U25.log", emit: log
    
    script:
    """
    export PATH=\$HOME/bioinformatics/UXM_deconv:\$PATH
    export PATH=\$HOME/bioinformatics/wgbs_tools:\$PATH
    
    echo "Running UXM deconvolution (U25) for ${sample_id}" | tee ${sample_id}_deconv_U25.log
    echo "Started: \$(date)" | tee -a ${sample_id}_deconv_U25.log
    
    uxm deconv \
        --atlas ${atlas} \
        --output ${sample_id}_deconv_U25.csv \
        ${pat} \
        2>&1 | tee -a ${sample_id}_deconv_U25.log
    
    echo "Completed: \$(date)" | tee -a ${sample_id}_deconv_U25.log
    """
}

/*
 * PROCESS 4: UXM Deconvolution - U250 Atlas
 */
process UXM_DECONV_U250 {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/deconvolution", mode: 'copy'
    
    cpus 4
    memory '16 GB'
    
    input:
    tuple val(sample_id), path(pat)
    path atlas
    
    output:
    tuple val(sample_id), path("${sample_id}_deconv_U250.csv"), emit: deconv
    path "${sample_id}_deconv_U250.log", emit: log
    
    script:
    """
    export PATH=\$HOME/bioinformatics/UXM_deconv:\$PATH
    export PATH=\$HOME/bioinformatics/wgbs_tools:\$PATH
    
    echo "Running UXM deconvolution (U250) for ${sample_id}" | tee ${sample_id}_deconv_U250.log
    echo "Started: \$(date)" | tee -a ${sample_id}_deconv_U250.log
    
    uxm deconv \
        --atlas ${atlas} \
        --output ${sample_id}_deconv_U250.csv \
        ${pat} \
        2>&1 | tee -a ${sample_id}_deconv_U250.log
    
    echo "Completed: \$(date)" | tee -a ${sample_id}_deconv_U250.log
    """
}

/*
 * PROCESS 5: Summarize and Compare Results
 */
process SUMMARIZE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    cpus 2
    memory '8 GB'
    
    input:
    path u25_files
    path u250_files
    
    output:
    path "all_samples_U25.csv", emit: u25_combined
    path "all_samples_U250.csv", emit: u250_combined
    path "summary_report.txt", emit: report
    
    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import glob
    
    # Combine U25 results
    u25_files = glob.glob("*_deconv_U25.csv")
    u25_dfs = []
    for file in u25_files:
        df = pd.read_csv(file)
        u25_dfs.append(df)
    
    if u25_dfs:
        u25_combined = pd.concat(u25_dfs, axis=1)
        u25_combined.to_csv("all_samples_U25.csv", index=False)
    
    # Combine U250 results
    u250_files = glob.glob("*_deconv_U250.csv")
    u250_dfs = []
    for file in u250_files:
        df = pd.read_csv(file)
        u250_dfs.append(df)
    
    if u250_dfs:
        u250_combined = pd.concat(u250_dfs, axis=1)
        u250_combined.to_csv("all_samples_U250.csv", index=False)
    
    # Create summary report
    with open("summary_report.txt", "w") as f:
        f.write("=== UXM Deconvolution Summary ===\\n\\n")
        f.write(f"Total samples analyzed: {len(u25_files)}\\n\\n")
        
        if u25_dfs:
            f.write("=== Top cell types across samples (U25 Atlas) ===\\n")
            # Get average fractions across samples
            numeric_cols = u25_combined.select_dtypes(include='number').columns
            if len(numeric_cols) > 0:
                means = u25_combined[numeric_cols].mean(axis=1)
                top_idx = means.nlargest(10).index
                f.write("Top 10 cell types by average fraction:\\n")
                for idx in top_idx:
                    cell_type = u25_combined.iloc[idx, 0]
                    avg = means.iloc[idx]
                    f.write(f"  {cell_type}: {avg:.4f}\\n")
        
        f.write("\\nAnalysis complete!\\n")
    
    print("Summary created successfully!")
    """
}

/*
 * WORKFLOW
 */
workflow {
    // Read sample sheet
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row.sample_id
            def bam_dir = row.bam_directory
            def bams = file("${bam_dir}/*.bam")
            tuple(sample_id, bams)
        }
    
    // Load atlases
    atlas_u25_ch = Channel.fromPath(params.atlas_u25)
    atlas_u250_ch = Channel.fromPath(params.atlas_u250)
    
    // Step 1: Merge BAMs
    MERGE_BAMS(samples_ch)
    
    // Step 2: Convert to PAT
    BAM_TO_PAT(MERGE_BAMS.out.merged)
    
    // Step 3: Deconvolve with both atlases
    UXM_DECONV_U25(BAM_TO_PAT.out.pat, atlas_u25_ch)
    UXM_DECONV_U250(BAM_TO_PAT.out.pat, atlas_u250_ch)
    
    // Step 4: Summarize results
    SUMMARIZE_RESULTS(
        UXM_DECONV_U25.out.deconv.map{ it[1] }.collect(),
        UXM_DECONV_U250.out.deconv.map{ it[1] }.collect()
    )
}

workflow.onComplete {
    log.info """\
        ========================================
        UXM Deconvolution Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Results   : ${params.outdir}
        ========================================
        """
        .stripIndent()
}
