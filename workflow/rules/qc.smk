
rule qc:
    input:
        r1 = "../fq/{sample}_R1_001.fastq.gz",
        r2 = "../fq/{sample}_R2_001.fastq.gz"
        #r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        #r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
    output:
        "results/qc/fastqc/{sample}_R1_001_fastqc.html",
        "results/qc/fastqc/{sample}_R1_001_fastqc.zip",
        "results/qc/fastqc/{sample}_R2_001_fastqc.html",
        "results/qc/fastqc/{sample}_R2_001_fastqc.zip"
    #log: "log/{sample}.qc.log.txt"
    #resources:
    #    mem = 1000,
    #    time = 300
    threads: 1
    message: """--- Quality check of raw data with FastQC before trimming."""
    shell: """
        fastqc -o results/qc/fastqc/ -f fastq {input.r1} &
        fastqc -o results/qc/fastqc/ -f fastq {input.r2}
        """





## RSEQC


rule rseqc_gtf2bed:
    input:
        "resources/genome.gtf",
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    #conda:
    #    "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/star/{sample}/aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}.log",
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}.log",
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/star/{sample}/aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}.log",
    #conda:
    #     "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/star/{sample}/aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/star/{sample}/aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    #conda:
    #    "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


#rule multiqc:
#    input:
#        expand(
#            "results/star/{unit.sample_name}/aligned.out.bam",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.junctionanno.junction.bed",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.junctionsat.junctionSaturation_plot.pdf",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.infer_experiment.txt",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.stats.txt",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.inner_distance_freq.inner_distance.txt",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.readdistribution.txt",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.readdup.DupRate_plot.pdf",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "results/qc/rseqc/{unit.sample_name}.readgc.GC_plot.pdf",
#            unit=units.itertuples(),
#        ),
#        expand(
#            "logs/rseqc/rseqc_junction_annotation/{unit.sample_name}.log",
#            unit=units.itertuples(),
#        ),
#    output:
#        "results/qc/multiqc_report.html",
#    log:
#        "logs/multiqc.log",
#    wrapper:
#        "v1.23.1/bio/multiqc"
