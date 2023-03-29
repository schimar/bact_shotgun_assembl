rule star_align:
    input:
        fq1 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R1_001.fastq.gz'])], 
        fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
        #fq1 = '/'.join([fqDir, '{sample}/{sample}_R1.fq.gz']), 
        #fq2 = '/'.join([fqDir, '{sample}/{sample}_R2.fq.gz']), 
        #fq1= "results/trimmed/{sample}_{unit}_R1.fastq.gz",
        #fq2= "results/trimmed/{sample}_{unit}_R2.fastq.gz",
        #unpack(get_fq),
        idx="resources/star_genome",
        gtf = "resources/genome.gtf",
    output:
        #directory("results/star/{sample}"),
        aln = "results/star/{sample}/aligned.out.bam",
        #log = "logs/star/{sample}/log.out",
        sj = "results/star/{sample}/SJ.out.tab",
        #aln="results/star/{sample}_{unit}/Aligned.out.bam",
        #sj="results/star/{sample}_{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}.log", #_{unit}.log",
    params:
        #index=lambda wc, input: input.index,
        extra="--chimSegmentMin 12 --chimOutType WithinBAM --readFilesSAMattrKeep All --quantMode GeneCounts --limitBAMsortRAM 24000000000 --outSAMtype BAM SortedByCoordinate --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30", # --sjdbGTFfile genome.gtf ",# --outSAMattrRGline ID:{sample} --sjdbGTFfile {} {}".format(
        #"resources/genome.gtf", config["params"]["star"]
        #),
    threads: 24
    #wrapper:
    #    "v1.23.4/bio/star/align"
    run:
      input_str_fq1 = ",".join(input.fq1)
      input_str_fq2 = ",".join(input.fq2)
      input_str = " ".join([input_str_fq1, input_str_fq2])
      shell(
        "STAR "
        " --runThreadN {threads}"
        " --genomeDir {input.idx}"
        " --readFilesIn {input_str}"
        " --readFilesCommand unpigz -c"
        " {params.extra}"
        #" --outTmpDir {tmpdir}/STARtmp"
        " --outFileNamePrefix results/star/{wildcards.sample}/"
        " --outStd BAM_SortedByCoordinate"
        " > {output.aln}"
        " {log}"
      )

# --quantMode GeneCounts 
# will count while mapping (NOTE: check if this doesn't need too much resources!!) 

rule samtools_sort:
    input:
        "results/star/{sample}_{unit}/aligned.out.bam",
    output:
        "results/samtools/{sample}_{unit}/alnd.sorted.bam",
    log:
        "logs/samtools/{sample}_{unit}.sort.log",
    params:
        extra="-m 4G",
    threads: 8
    #wrapper:
    #    "v1.23.4/bio/samtools/sort"
    shell:
      """
      samtools sort -@ {threads} -o {output} {input} > {log} 2>&1
      """


rule samtools_merge:
    input:
        ["results/samtools/{sample}_L001/alnd.sorted.bam", "results/samtools/{sample}_L002/alnd.sorted.bam"],
    output:
        "results/samtools/{sample}.mrgd.bam",
    log:
        "logs/samtools/{sample}_merged.log",
    params:
        extra="",  # optional additional parameters as string
    threads: 8
    #wrapper:
    #    "v1.23.4/bio/samtools/merge"
    shell:
      """
      samtools merge -@ {threads} -o {output} {input} > {log} 2>&1
      """


rule samtools_flagstat:
    input:
        "results/star/{sample}/aligned.out.bam",
    output:
        "results/samtools/flagstat/{sample}.bam.flagstat",
    log:
        "logs/samtools/flagstat/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.23.5/bio/samtools/flagstat"
        


rule rRNA_bwa_mem:
    input:
        fq1 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R1_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R1_001.fastq.gz'])], 
        fq2 = ['/'.join([fqDir, 'Fastq/{sample}_L001_R2_001.fastq.gz']), '/'.join([fqDir, 'Fastq/{sample}_L002_R2_001.fastq.gz'])], 
        idx = multiext("resources/genome.fa", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        "results/bwa_mem2/{sample}.rRNA.bam",
    log:
        "logs/bwa_mem2/{sample}.log",
    params:
        extra = r"-R '@RG\tID:{sample}'",
        sort = "samtools", # "none", "samtools" or "picard"
        sort_order = "coordinate",  # "coordinate" (default) or "queryname"
        sort_extra = "",
    threads: 12
    wrapper:
        "v1.23.5/bwa-mem2/mem" 


