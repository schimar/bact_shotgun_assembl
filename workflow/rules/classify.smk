rule kraken:
    input:
        reads = "results/spades/{sample}/contigs.fasta",
        #reads = lambda wildcards: sample_reads[wildcards.samp],
    output:
        krak = "results/kraken/{sample}.krak",
        krak_report = "results/kraken/{sample}.krak.report",
    params:
        db = config['kraken_db'],
        paired_string = '',  #'--paired',
        confidence_threshold = 0.0
    threads: 28
    resources:
        mem = 200,
        #time=6
    #singularity: "docker://quay.io/biocontainers/kraken2:2.1.2--pl5262h7d875b9_0"
    # singularity: "docker://quay.io/biocontainers/kraken2:2.0.9beta--pl526hc9558a2_0"
    shell: 
        """
        time kraken2 --db {params.db} --threads {threads} --output {output.krak} \
        --report {output.krak_report} {params.paired_string} {input.reads} \
        --confidence {params.confidence_threshold} --use-names 
        """
# --fastq-input --gzip-compressed

rule bracken:
    input:
        krak_report = "results/kraken/{sample}.krak.report",
        krak = "results/kraken/{sample}.krak",
    output:
        #"results/bracken/{sample}.krak_bracken_species.report",
        "results/bracken/{sample}.krak.report.bracken",
    params:
        db = config['kraken_db'],
        readlen = 150,   #config['read_length'],
        level = "S",    #config['taxonomic_level'],
        threshold = 10,
        outspec = "results/bracken/{sample}.krak.report.bracken",
    threads: 1
    resources:
        mem = 64,
        #time = 1
    shell: 
        """
        bracken -d {params.db} -i {input.krak_report} -o {params.outspec} -r {params.readlen} \
        -l {params.level} -t {params.threshold}
        """

rule krona:
    input: rules.kraken.output.krak_report,
    output: "results/krona/{sample}.html",
    params:
        db = config['kraken_db'],
    shell: 
        """
        ktImportTaxonomy -q 2 -t 3 {input} -o {output}
                """
#ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i {input} -o {output} -tax {params.db}

        
#-tax $(which kraken2 | sed 's/envs\/classification2.*$//g')/envs/classification2/bin/taxonomy



