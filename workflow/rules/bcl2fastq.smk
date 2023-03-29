rule bcl2fastq:
    input:
        rundir = "/mnt/illumina/230209_A01272_0035_BHTVFGDRX2/",
    output:
        expand("results/bcl2fq/{runid}/{sample}_{lane}_{read}_001.fastq.gz", sample = list(samples_dict.keys()), lane = ['L001', 'L002'], read = ['R1', 'R2'], runid = config['runID']),
    params:
        outdir = "results/bcl2fq/",
    #log:
    #    "logs/bcl2fastq.log",
    threads: 24
    shell:
        """
        nohup bcl2fastq --runfolder-dir {input.rundir} --output-dir {params.outdir} -p 18 -r 3 -w 3 > logs/bcl2fastq.log 2>&1
        """
# expand("../fq/{{runid}}_S0_L001_{readid}_001.fastq.gz", readid=config['readids'])
# I'll probably have to use the expand in the output, so I can name the files in rule all! 
    
#   "{sample}_{read}_001.fastq.gz"
