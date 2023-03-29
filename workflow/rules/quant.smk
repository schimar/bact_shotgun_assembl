
rule salmon_decoys:
    input:
        transcriptome="resources/gencode.v43.transcripts.fa",
        genome="resources/genome.fa",
    output:
        gentrome="resources/gentrome.fa",
        decoys="resources/decoys.txt",
    threads: 2
    log:
        "logs/salmon/decoys.log",
    wrapper:
        "v1.23.5/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences="resources/gentrome.fa",
    output:
      #multiext(
        directory("resources/gentrome_index/"),
    log:
        "logs/salmon/index.log",
    threads: 6
    params:
        # optional parameters
        extra="",
    wrapper:
        "v1.23.5/bio/salmon/index"

  
rule salmon_quant:
    input:
      fq1 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz'])], 
      fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L002_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
      #r1="results/trimmed/{sample}_R1.fastq.gz",
      #r2="results/trimmed/{sample}_R2.fastq.gz",
      index="resources/gentrome_index",
    output:
      quant="results/salmon/{sample}/quant.sf",
      lib="results/salmon/{sample}/lib_format_counts.json",
    log:
        "logs/salmon/{sample}/quant.log",
    threads: 6
    wrapper:
        "v1.23.5/bio/salmon/quant"



rule htseq:
        input:
            bam = 'results/star/{sample}/aligned.out.bam',
            gtf = 'resources/genome.gtf'
        output:
            union = 'logs/htseq/{sample}_HTSeq_union_gtf_no_gene_ID.log',
            csv = 'results/htseq/{sample}_HTSeq.csv'
        log:
            'logs/htseq/{sample}.log',
        threads: 1
        shell:
            'htseq-count -m union -s no -t gene -i gene_id -r pos -f bam {input.bam} {input.gtf} &> {output.union} && '
            'grep ENS {output.union} | sed "s/gene://g" > {output.csv}'



