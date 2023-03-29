rule arriba:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        #bam="results/star/{sample}/aligned.out.bam",
        genome="resources/genome.fa",
        annotation="resources/genome.gtf",
        # optional: # A custom tsv containing identified artifacts, such as read-through fusions of neighbouring genes.
        #blacklist="resources/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
        #known_fusions="resources/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
        # default blacklists are selected via blacklist parameter
        # see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
        #custom_blacklist=[],
    output:
        fusions="results/arriba/{sample}/fusions.tsv",
        discarded="results/arriba/{sample}/fusions.discarded.tsv",
        #done="results/arriba/{sample}/arriba.done",
    params:
        # required if blacklist or known_fusions is set
        genome_build="GRCh38",
        default_blacklist=True,
        default_known_fusions=True,
        extra="",   #alignIntronMax",
    log:
        "logs/arriba/{sample}.log",
    threads: 1
    wrapper:
        "v1.23.4/bio/arriba"


rule catL_fqs:
    input: 
        fq1 = '/'.join([fqDir, 'Fastq/{sample}_L001_{read}_001.fastq.gz']), 
        fq2 = '/'.join([fqDir, 'Fastq/{sample}_L002_{read}_001.fastq.gz']), 
        #fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
    output:
        cat_fq = "results/fq_cats/{sample}_{read}_catL.fq.gz",
    log:
        "logs/catL/{sample}_{read}.log",
    threads: 10
    shell:
      """
        zcat {input.fq1} {input.fq2} | pigz -p {threads} -c > {output.cat_fq} 2> {log}
      """



# cat fq.gz files before? 
rule starseqr:
    input: 
      fq1="results/trimmed/{sample}_{unit}_R1.fastq.gz",
      fq2="results/trimmed/{sample}_{unit}_R2.fastq.gz",
      genome="resources/genome.fa",
      genomeDir="resources/star_genome",
      gtf="resources/genome.gtf",
    output:
      directory("results/starseqr/{sample}/") 
    shell:
      """
      starseqr.py -1 {fq1} -2 {fq2} -m 1 -p  -t 12 -i {genomeDir} -g {gtf} -r {genome} -vv
      """



rule star_fusion:
    input:
      fq
    output:

