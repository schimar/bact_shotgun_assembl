rule vcftls_relate2:
  input:
    vcf = 'vars/ta{vartype}Bypop.vcf',
  output:
    touch('vars/ta{vartype}/stats/vcftools/king.done') 
  message:
    """--- Calculating relatedness2 (Manichaikul) with vcftools ---"""
  shell:
    """
    vcftools --vcf {input.vcf} --relatedness2 --out vars/ta{wildcards.vartype}/stats/vcftools/king
    """

rule relatePlot:
  input:
    zarr = 'vars/ta{type}.zarr/',
    rel2 = 'vars/ta{type}/stats/vcftools/king.relatedness2'
  output:
    touch('vars/ta{type}/figs/vcftools/plots.done')
  message:
    """--- Plotting relatedness plots ---"""
  shell:
    """
    script/relatePlot.py {input.zarr} {input.rel2}
    """


#### Tajima's D  (w/ 1 kbp bin size)
#vcftools --gzvcf ../g12f3.vcf.gz --TajimaD 1000 --keep wk.txt --out tajD/wk_tD
#vcftools --gzvcf ../g12f3.vcf.gz --TajimaD 1000 --keep wp.txt --out tajD/wp_tD

# site-pi, Tajima's D, Fis, Fst, LROH, heterozygosity? 
