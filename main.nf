nextflow.enable.dsl=2

// Define parameters
params.chromosomes = (1..22).toList()
params.vcf_dir     = "/hpctmp/YOUR_ID/5004/data/vcf"
params.pheno       = "${projectDir}/data/eas_eur_phenotype.txt"
params.outdir      = "${projectDir}/results"

// Define processes
process run_qc {
    tag "chr${chr_idx}"
    label 'scaling_target'
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1'

    input:
    val chr_idx

    output:
    tuple val(chr_idx), path("chr${chr_idx}_qc.*")

    script:
    """
    plink \\
        --vcf ${params.vcf_dir}/ALL.chr${chr_idx}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \\
        --keep ${params.pheno} \\
        --make-bed \\
        --maf 0.05 \\
        --geno 0.05 \\
        --hwe 1e-6 \\
        --out chr${chr_idx}_qc
    """
}

process run_assoc {
    tag "chr${chr_idx}"
    label 'scaling_target'
    container 'docker://quay.io/biocontainers/plink:1.90b6.21--h779adbc_1'

    input:
    tuple val(chr_idx), path(qc_files)

    output:
    path "chr${chr_idx}_results.assoc"

    script:
    """
    plink \\
        --bfile chr${chr_idx}_qc \\
        --pheno ${params.pheno} \\
        --assoc \\
        --allow-no-sex \\
        --out chr${chr_idx}_results
    """
}

process merge_results {
    publishDir "${params.outdir}", mode: 'copy'
    container 'docker://quay.io/jupyter/scipy-notebook:latest'

    input:
    path assoc_files

    output:
    path "all_chromosomes.assoc"

    script:
    """
    python3 -c "
import pandas as pd, glob
files = sorted(glob.glob('*_results.assoc'))
df = pd.concat([pd.read_csv(f, sep=r'\s+', engine='python') for f in files])
df.to_csv('all_chromosomes.assoc', sep='\\t', index=False)
    "
    """
}

process plot_manhattan {
    publishDir "${params.outdir}", mode: 'copy'
    container 'docker://quay.io/jupyter/scipy-notebook:latest'

    input:
    path merged_file

    output:
    path "manhattan_all.png"

    script:
    """
    python3 ${projectDir}/scripts/plot_manhattan.py \\
        --input ${merged_file} \\
        --out manhattan_all.png
    """
}

// Define the workflow
workflow {
    ch_chromosomes = Channel.fromList(params.chromosomes)

    ch_qc = run_qc(ch_chromosomes)

    ch_assoc = run_assoc(ch_qc)

    ch_merged = merge_results(ch_assoc.collect())

    plot_manhattan(ch_merged)
}