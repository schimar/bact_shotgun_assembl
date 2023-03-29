#rule trim:
#    input:
#        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
#        r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
#        adapters = config["adapters"]
#    output:
#        r1trmd = "trm/{sample}_R1.fq.gz",
#        r2trmd = "trm/{sample}_R2.fq.gz"
#    threads: 2
#    message: """--- Quality trimming of fastq files before mapping."""
#    shell: 
#        """
#        bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1trmd} out2={output.r2trmd} trimq=6 qtrim=r hdist=1 bhist=trm/hist/{wildcards.sample}.bhist qhist=trm/hist/{wildcards.sample}.qhist lhist=trm/hist/{wildcards.sample}.lhist tpe tbo 
#        """


rule bbduk_pe:
    input:
        sample =["../fq/{sample}_R1_001.fastq.gz", "../fq/{sample}_R2_001.fastq.gz"],
        adapters ="resources/adapters.fa",   # from bbmap/resources/
    output:
        trimmed =["results/trimmed/{sample}.1.fq.gz", "results/trimmed/{sample}.2.fq.gz"],
        singleton ="results/trimmed/{sample}.single.fq.gz",
        discarded ="results/trimmed/{sample}.discarded.fq.gz",
        stats ="results/trimmed/{sample}.stats.txt",
    wildcard_constraints:
        sample = '[A-Za-z0-9\_\-]+L00[0-9A-Z\-\_.a-z]+'
    log:
        "logs/bbduk/{sample}.log"
    params:
        extra = lambda w, input: "ref={},adapters,artifacts ktrim=r k=23 mink=11 hdist=1 tpe tbo trimpolygright=10 minlen=25 maxns=30 entropy=0.5 entropywindow=50 entropyk=5".format(input.adapters),
    threads: 7
    wrapper:
        "v1.23.5/bio/bbtools/bbduk"

rule cat_fq1s:
    #input:
    #    "results/bbnorm/bbdone"
    output:
        "results/fqcat/{sample}.1.fq.gz",
    threads: 10
    shell:
        """
        zcat results/bbnorm/{wildcards.sample}_L00*_1nrm.fq.gz | pigz -c > {output} & 
        """

rule cat_fq2s:
    #input:
    #    "results/bbnorm/bbdone"   
    output:
        "results/fqcat/{sample}.2.fq.gz",
    threads: 10
    shell:
        """
        zcat results/bbnorm/{wildcards.sample}_L00*_2nrm.fq.gz | pigz -c > {output}
        """



