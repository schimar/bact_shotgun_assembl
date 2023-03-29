rule qc:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
    output:
        "raw/qc/fastqc/{sample}_R1_001_fastqc.html",
        "raw/qc/fastqc/{sample}_R1_001_fastqc.zip",
        "raw/qc/fastqc/{sample}_R2_001_fastqc.html",
        "raw/qc/fastqc/{sample}_R2_001_fastqc.zip"
    #log: "log/{sample}.qc.log.txt"
    #resources:
    #    mem = 1000,
    #    time = 300
    threads: 1
    message: """--- Quality check of raw data with FastQC before trimming."""
    shell: """
        fastqc -o raw/qc/fastqc/ -f fastq {input.r1} &
        fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
        """


rule trim:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
        adapters = config["adapters"]
    output:
        r1trmd = "trm/{sample}_R1.fq.gz",
        r2trmd = "trm/{sample}_R2.fq.gz"
    threads: 2
    message: """--- Quality trimming of fastq files before mapping."""
    shell: 
        """
        bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1trmd} out2={output.r2trmd} trimq=6 qtrim=r hdist=1 bhist=trm/hist/{wildcards.sample}.bhist qhist=trm/hist/{wildcards.sample}.qhist lhist=trm/hist/{wildcards.sample}.lhist tpe tbo 
        """

rule refIndex:
	input:
		ref = config['ref']
	output:
		'ref/genome/1/summary.txt'
	shell:
		"""
		bbmap.sh ref={input.ref}
		"""


rule map:
	input:
		tr1 = lambda wildcards: getTrmHome(wildcards.sample)[0],
		tr2 = lambda wildcards: getTrmHome(wildcards.sample)[1],
		ref = config["ref"]
	output:
		bam = "bbmap/{sample}.bam"
	threads: 24
	message: """--- Mapping reads to reference genome ---"""
	shell:
		"""
		bbmap.sh -Xmx50g t={threads} ref={input.ref} in1={input.tr1} in2={input.tr2} out=stdout.sam minid=0.85 rgid={wildcards.sample} rglb=igaDNA rgsm={wildcards.sample} rgpl=ILLUMINA overwrite=f unpigz=t | samtools view -F 4 -Shu - | samtools sort - -o {output.bam}
		"""

rule bamIndex:
	input: 
		aln = 'bbmap/{sample}.bam'	
	output:
		idx = "bbmap/{sample}.bam.bai"
	threads: 2
	message: """--- Indexing with samtools ---"""
	shell:
		"""
		samtools index {input.aln} {output.idx}
		"""





