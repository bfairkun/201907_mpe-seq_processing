genomeDir = config["genome_dir"][:-1]
rule STAR_make_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["genome_gtf"]
    output:
        config["genome_dir"] + "chrLength.txt",
    log:
        "logs/STAR/MakingIndex.log"
    params:
        genomeDir = genomeDir
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule umi_tools_extract:
    input:
        fastq=GetRNASeqFastq
    output:
        R1 = output_dir + "umi_extract_fastq/{sample}/R1.fastq.gz",
        R2 = output_dir + "umi_extract_fastq/{sample}/R2.fastq.gz"
    log:
        output_dir + "logs/umi_tools_extract/{sample}.log"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        conda activate umi_tools
        umi_tools extract -I {input[0]} --bc-pattern=NNNNNNN --read2-in={input[1]} --stdout={output.R1} --read2-out={output.R2} &> {log}
        """

rule STAR_alignment:
    input:
        index = config["genome_dir"] + "chrLength.txt",
        R1 = output_dir + "umi_extract_fastq/{sample}/R1.fastq.gz",
        R2 = output_dir + "umi_extract_fastq/{sample}/R2.fastq.gz"
    log:
        output_dir + "logs/STAR/{sample}.log"
    threads: 8
    output:
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir {genomeDir} --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --outWigStrand Stranded --outWigType wiggle --alignEndsType EndToEnd --readFilesCommand zcat --outFileNamePrefix {output_dir}Alignments/{wildcards.sample}/ --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --alignSplicedMateMapLmin 18 &> {log}
        # run ulimit -v 30953070230
        # STAR --runThreadN {threads} --genomeDir /project2/yangili1/RNAseq_pipeline/index/GRCh37/STAR_hg19/ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --outWigStrand Stranded --outWigType wiggle --alignEndsType EndToEnd --readFilesCommand zcat --outFileNamePrefix {output_dir}Alignments/{wildcards.sample}/ --clip3pAdapterSeq CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --alignSplicedMateMapLmin 18 --outFilterIntronMotifs RemoveNoncanonicalUnannotated &> {log}
        """

rule index_RNA_seq_bams:
    input:
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        output_dir + "logs/index_RNA_seq_bams/{sample}.log"
    shell:
        "samtools index {input} &> {log}"

rule umi_tools_dedup_saturation_analysis:
    input:
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        output_dir + "Dedup_UMISaturationAnalysis/{sample}_{subset_frac}.bam"
    log:
        output_dir + "logs/umi_tools_dedup_saturation_analysis/{sample}.{subset_frac}.log"
    conda:
        "../envs/umi_tools.yaml"
    shell:
        """
        umi_tools dedup -I {input} --paired -S {output} --spliced-is-unique --method unique --unmapped-reads use --subset {wildcards.subset_frac} &> {log}
        """

rule umi_tools_dedup:
    input:
        output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bam = output_dir + "Dedup/{sample}.bam",
        bai = output_dir + "Dedup/{sample}.bam.bai",
    log:
        output_dir + "Dedup/{sample}.log"
    conda:
        "../envs/umi_tools.yaml"
    params:
        subsample = 0.16
    shell:
        """
        umi_tools dedup -I {input} --paired -S {output.bam} --spliced-is-unique --method unique --unmapped-reads use --subset {params.subsample} &> {log}
        samtools index {output.bam}
        """

rule remove_unextended_primer:
    input:
        bam = output_dir + "Dedup/{sample}.bam",
        bai = output_dir + "Dedup/{sample}.bam.bai",
    output:
        bam = output_dir + "Dedup/{sample}.extensions.bam",
        bai = output_dir + "Dedup/{sample}.extensions.bam.bai",
    shell:
        """
        cat <(samtools view -H {input.bam}) <(samtools view {input.bam} | perl -F'\\t' -lane 'BEGIN {use List::Util qw(sum);} @mymatches = $F[5] =~ m/(\d+)M/g; if( (abs($F[8]) >= 30) || (sum(@mymatches)>30)  ) {print}') | samtools view -bS - > {output.bam}
        samtools index {output.bam}
        """

rule CountJunctions_Orig:
    input:
        bam = output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bai = output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        junc = output_dir + "Alignments/{sample}/Aligned.sortedByCoord.out.bam.junc",
    log:
        output_dir + "logs/CountJunctions_Orig/{sample}.log"
    shell:
        """
        bam2junc.sh {input.bam} {output.junc}
        """

rule CountJunctions_Dedup:
    input:
        bam = output_dir + "Dedup/{sample}.bam",
        bai = output_dir + "Dedup/{sample}.bam.bai",
    output:
        junc = output_dir + "Dedup/{sample}.bam.junc",
    log:
        output_dir + "logs/CountJunctions_Dedup/{sample}.log"
    shell:
        """
        bam2junc.sh {input.bam} {output.junc}
        """

#Get annotated splice sites (transcripts minus exons)
# bedtools subtract -s -a <(awk -F'\t' '$3=="transcript"' /project2/gilad/bjf79/genomes/GRCh38_Ensembl/Annotations/Homo_sapiens.GRCh38.94.chr.gtf | bedtools sort -i - -g /project2/gilad/bjf79/genomes/GRCh38_Ensembl/Sequence/Homo_sapiens.GRCh38.dna_sm.chromosome.all.fa.fai) -b <(awk -F'\t' '$3=="exon"' /project2/gilad/bjf79/genomes/GRCh38_Ensembl/Annotations/Homo_sapiens.GRCh38.94.chr.gtf | bedtools sort -i - -g /project2/gilad/bjf79/genomes/GRCh38_Ensembl/Sequence/Homo_sapiens.GRCh38.dna_sm.chromosome.all.fa.fai) > Output/scratch/Introns.bed

#Get 3ss
# cat Introns.bed | awk -F'\t' '$7=="+" {print $1"_"$5+1"_+"} $7=="-" {print $1"_"$4"_-"}' | sort | uniq > Annotated3ss.txt

# Filter targets by annotated3ss

#convert to junc with leafcutter script, strand is opposite actual strand

#Convert to count by 3ss, where duplicates are allowed
