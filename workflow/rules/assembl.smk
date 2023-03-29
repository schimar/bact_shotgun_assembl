#ruleorder: bbnorm > cat_fq1s > cat_fq2s

rule bbnorm:
    input:
        #fq1 = "results/trimmed/{sample}.1.fq.gz",
        #fq2 = "results/trimmed/{sample}.2.fq.gz",
        fq1 = "/home/schimar/bio/ref/bact/Streptomyces_netropsis/SRR10053378/{sample}_1.fq.gz",
        fq2 = "/home/schimar/bio/ref/bact/Streptomyces_netropsis/SRR10053378/{sample}_2.fq.gz",
    output:
        nrm1 = "results/bbnorm/{sample}_1nrm.fq.gz",
        nrm2 = "results/bbnorm/{sample}_2nrm.fq.gz",
        #nrm1 = "/home/schimar/bio/ref/bact/Streptomyces_netropsis/SRR10053378/{sample}_1nrm.fq.gz",
        #nrm2 = "/home/schimar/bio/ref/bact/Streptomyces_netropsis/SRR10053378/{sample}_2nrm.fq.gz",  
    wildcard_constraints:
        sample = '[A-Za-z0-9\_\-]+'   #L00[0-9A-Z\-\_.a-z]+'
    params:
        targetCov = 40,
        minDepth = 2,
    threads: 16
    log:
        "logs/bbnorm/{sample}.log",
    shell:
        """
        bbnorm.sh -Xmx30g -threads={threads} in={input.fq1} in2={input.fq2} target={params.targetCov} mindepth={params.minDepth} out={output.nrm1} out2={output.nrm2} > {log} 2>&1
        """

rule spades:
    input:
        #fq1 = "results/fqcat/{sample}.1.fq.gz",
        #fq2 = "results/fqcat/{sample}.2.fq.gz",
        fq1 = "results/bbnorm/{sample}_1nrm.fq.gz",
        fq2 = "results/bbnorm/{sample}_2nrm.fq.gz",
    output:
        pth = directory("results/spades/{sample}/"),
        #done = touch("results/spades/{sample}/spades.done"),
        fa = "results/spades/{sample}/contigs.fasta",
    log:
        "logs/spades/{sample}.log",
    wildcard_constraints:
        sample = '[A-Za-z0-9\_\-]+'   #L00[0-9A-Z\-\_.a-z]+'
    threads: 32
    shell:
        """
        spades.py -t {threads} --isolate -1 {input.fq1} -2 {input.fq2} -o {output.pth} > {log} 2>&1
        """
#-k 21,33,55,77 --careful
#          #(' --trusted-contigs ' + 'data/input/' + config['reference_genome'] + '/' + config['genome_file_identifier']))
## phred-offset 33
# if using 2x150bp illumina reads, see https://cab.spbu.ru/files/release3.14.0/manual.html#sec3.4
# -k 21,33,55,77

rule quast:
    input:
        fasta="results/spades/{sample}/contigs.fasta",
#        #ref="genome.fasta",
#        #gff="annotations.gff",
        #pe1 = "results/fqcat/{sample}.1.fq.gz",
        #pe2 = "results/fqcat/{sample}.2.fq.gz",
        pe1 = "results/bbnorm/{sample}_1nrm.fq.gz",
        pe2 = "results/bbnorm/{sample}_2nrm.fq.gz",#        #pe12="reads.fastq",
#        #mp1="matereads_R1.fastq",
#        #mp2="matereads_R2.fastq",
#        #mp12="matereads.fastq",
#        #single="single.fastq",
#        ##pacbio="pacbio.fas",
#        ##nanopore="nanopore.fastq",
#        ##ref_bam="ref.bam",
#        ##ref_sam="ref.sam",
#        #bam=["s1.bam","s2.bam"],
#        #sam=["s1.sam","s2.sam"],
#        #sv_bedpe="sv.bed",
    output:
        multiext("results/quast/{sample}/report.", "html", "tex", "txt", "pdf", "tsv"),
        multiext("results/quast/{sample}/transposed_report.", "tex", "txt", "tsv"),
        #multiext(
        #    "results/quast/{sample}/basic_stats/",
        #    "cumulative_plot.pdf",
        #    "GC_content_plot.pdf",
        #    "gc.icarus.txt",
        #    "genome_GC_content_plot.pdf",
        #    "NGx_plot.pdf",
        #    "Nx_plot.pdf",
        #),
        #multiext(
        #    "results/quast/{sample}/contigs_reports/",
        #    "all_alignments_genome.tsv",
        #    "contigs_report_genome.mis_contigs.info",
        #    "contigs_report_genome.stderr",
        #    "contigs_report_genome.stdout",
        #),
        #"results/quast/{sample}/contigs_reports/minimap_output/genome.coords_tmp",
        "results/quast/{sample}/icarus.html",
        "results/quast/{sample}/icarus_viewers/contig_size_viewer.html",
        "results/quast/{sample}/quast.log",
    log:
        "logs/quast/{sample}.quast.log",
    params:
        extra="--min-contig 5 --min-identity 95.0",
    wrapper:
        "v1.23.5/bio/quast"
# quast.py Hcov.contigs.fasta -1 ~/bio/ref/bact/Streptomyces_netropsis/SRR10053378/SRR10053378_1nrm.fq.gz -2 ~/bio/ref/bact/Streptomyces_netropsis/SRR10053378/SRR10053378_2nrm.fq.gz -o quast_Hcov

rule filter_fasta:
    input:
        fa = "results/spades/{sample}/contigs.fasta",
    output:
        fa = "results/spades/{sample}/contigs.hcov.fa",
    shell:
        """
        for file in {input.fa}; do
          grep -F '>' $file | sed -e 's/_/ /g' | sort -nrk 6 | awk '$6>=4.0 && $4>=200 {{print $0}}' | sed -e 's/ /_/g' | sed -e 's/>//g' > $file.txt
          echo sequences to keep
          wc -l $file.txt 
          echo running fastagrep.pl
          scripts/fastagrep.pl -f $file.txt $file > {output.fa}   #$filebase.hcov.fa
          echo sequences kept
          grep -c '>' {output.fa}
        done
        """


rule prokka:
    input:
        fa = "results/spades/{sample}/contigs.fasta",
    output:
        pth = directory("results/prokka/{sample}"),
        #done = touch("results/prokka/{sample}/prokka.done"),
        gff = "results/prokka/{sample}/prok.gff",
    wildcard_constraints:
        sample = '[A-Za-z0-9\_\-]+'   #L00[0-9A-Z\-\_.a-z]+'
    log:
        "logs/prokka/{sample}.log",
    shell:
        """
        prokka {input.fa} --force --outdir {output.pth} -prefix prok > {log} 2>&1
        """

rule clean_gff:
    input:
        "results/prokka/{sample}/prok.gff",
    output:
        "results/prokka/{sample}/prok_noSeq.gff"
    log:
        "logs/prokka/{sample}_clngff.log",
    shell:
        """
        sed '/^##FASTA/Q' {input} > {output}
        """


