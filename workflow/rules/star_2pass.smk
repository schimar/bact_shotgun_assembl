# adapted from https://evodify.com/rna-seq-star-snakemake/

#'/'.join([fqDir, '{sample}_L001_R1_001.fastq.gz'])

rule star_p1:
        input:
            fq1 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R1_001.fastq.gz'])], 
            fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
            # paired end reads needs to be ordered so each item in the two lists match
            #fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
            idx = "resources/star_genome",
            #gtf = "resources/genome.gtf",
        output:
            aln = "results/star/{sample}_p1/aligned.out.bam",
            log = "logs/star/{sample}_p1/log.out",
            sj = "results/star/{sample}_p1/sj.out.tab",
        log:
            "logs/star/{sample}_p1.log",
        params:
            #index=lambda wc, input: input.index,
            extra="--chimSegmentMin 20 --chimOutType WithinBAM --readFilesSAMattrKeep All --limitBAMsortRAM 24000000000 --outSAMtype BAM Unsorted",
        threads: 20
        wrapper:
            "v1.23.5/bio/star/align"
# --sjdbGTFfile {input.gtf} 

            
rule sjdir:
        output:
            directory('results/star/sj'),
        threads: 1
        shell:
            'mkdir {output}'

# We filter all splice junctions with >= 3 supporting reads 
filt_cmd = r"""{if ($7 >= 3) print $0}"""

rule filter:
        input:
            'results/star/{sample}_p1/sj.out.tab',
            #directory('sj')
        output:
            'results/star/sj/{sample}_p1/sj.filtered.tab'
        log:
            'logs/star/{sample}_filter.log',
        threads: 1
        shell:
            """
            awk {filt_cmd:q} {input} > {output} 2> {log} 
            """
# {input}.filtered 2> {log} && 'mv {input}.filtered {output}'


#awk '{ { if (\$7 >= 3) print \$0 } }' {{input[0]}} > {input[0]}.filtered && 

rule star_p2:
    input:
        fq1 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R1_001.fastq.gz'])], 
        fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
        sjfltrd = 'results/star/sj/{sample}_p1/sj.filtered.tab',
        idx = "resources/star_genome",
        gtf = "resources/genome.gtf",
    output:
        aln = "results/star/{sample}_p2/aligned.sortedByCoord.out.bam",
        log = "logs/star/{sample}_p2/log.out",
        sj = "results/star/{sample}_p2/sj.out.tab",
    log:
        "logs/star/{sample}_p2.log",
    params:
        #index=lambda wc, input: input.index,
        ids = '{sample}',
        extra="--chimSegmentMin 20 --chimOutType WithinBAM --readFilesSAMattrKeep All --limitBAMsortRAM 24000000000 --sjdbFileChrStartEnd results/star/sj/{sample}_p1/sj.filtered.tab --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:{sample}",
    threads: 20
    wrapper:
        "v1.23.5/bio/star/align"

# --sjdbGTFfile {input.gtf} 


rule samtools_flagstat:
    input:
        "results/star/{sample}_p2/aligned.sortedByCoord.out.bam",
    output:
        "results/samtools/flagstat/{sample}_p2/alnd.srtd.bam.flagstat",
    log:
        "logs/samtools/flagstat/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.23.4/bio/samtools/flagstat"
        



