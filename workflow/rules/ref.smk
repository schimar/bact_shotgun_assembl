rule get_genome:
    output:
        "resources/genome.fa",
    log:
        "logs/get_genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fa",
    output:
        "resources/genome.fa.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/genome.fa",
    output:
        multiext("resources/genome.fa", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        "logs/bwa-mem2_index/genome.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v1.23.5/bio/bwa-mem2/index"


rule star_index:
    input:
        fasta="resources/genome.fa",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 8
    params:
      extra=lambda wc, input:"--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "v1.23.4/bio/star/index"  
        #"0.77.0/bio/star/index"
