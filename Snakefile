# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
shell.executable("/bin/bash")
include: "rules/common.smk"
configfile: "config.yaml"

AllRNASeq = expand (output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam", sample=samples.index ) + expand (output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=samples.index )

subset = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1]
AllSaturationAnalysis = expand(output_dir + "Dedup_UMISaturationAnalysis/{sample}_{subset_frac}.bam", sample=samples.index, subset_frac=subset)


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        config["genome_dir"] + "chrLength.txt",
        expand(output_dir + "umi_extract_fastq/{sample}/R1.fastq.gz", sample=samples.index),
        # AllRNASeq,
        AllSaturationAnalysis,
        expand (output_dir + "Dedup/{sample}.bam.junc", sample=samples.index),
        expand (output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.junc", sample=samples.index)

include: "rules/other.smk"
include: "rules/STAR_Alignment.smk"
