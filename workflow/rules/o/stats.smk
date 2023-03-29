rule vcf2zarr:
  input:
    vcf = 'vars/ta{vartype}Bypop.vcf'
  output:
    'vars/ta{vartype}.zarr/.zgroup'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
      python script/vcf2zarr.py {input.vcf}
    """

##echo '{input.vcf}.zarr'


rule alStats:
  input:
    zarr = 'vars/ta{vartype}.zarr/'
  output:
    touch("vars/ta{vartype}/al.done")
  message:
    """--- Creating scikit-allel statistics ---"""
  shell:
    """
    python script/al.py {input.zarr}
    """

rule gemma_mg:
  input:
    vcf = 'vars/taSubInDel.{sets}.vcf'
    #vcf = 'vars/ta{vartype}Bypop.vcf'
  output:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg'
    #mg = 'vars/ta{vartype}/stats/gemma/{vartype}.mg'
  message:
    """--- Converting vcf to bimbam format ---"""
  shell:
    """
    python script/vcf2mg.py {input.vcf} > {output.mg}
    """


#Which of the two relatedness matrix to choose will largely depend on the underlying genetic architecture of the given trait. Specifically, if SNPs with lower minor allele frequency tend to have larger effects (which is inversely proportional to its genotype variance), then the standardized genotype matrix is preferred. If the SNP effect size does not depend on its minor allele frequency, then the centered genotype matrix is preferred. In our previous experience based on a limited examples, we typically find the centered genotype matrix provides better control for population structure in lower organisms, and the two matrices seem to perform similarly in humans.
rule gemma_relmat_stdzd:
  input:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg',
    pheno = 'vars/taSubInDel/stats/gemma/{trait}.pheno'
    #mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg',
    #pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    relM = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.stdzd.relmat'
    #relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.stdzd.relmat'
  message:
    """--- Calculating standardized relatedness matrix with gemma ---"""
  shell:
    """
    gemma -g {input.mg} -p {input.pheno} -gk 2 -o {output.relM}
    """

rule gemma_relmat_ctrd:
  # NOTE: see **manual 4.4.2** on details regarding centered vs standardized rel.matrix
  input:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg',
    pheno = 'vars/taSubInDel/stats/gemma/{trait}.pheno'
    #mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg',
    #pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    touch('vars/taSubInDel/stats/gemma/relmat/ctrd.done'),
    relM = 'taSubInDel'
        #relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.ctrd.relmat'
  message:
    """--- Calculating centered relatedness matrix with gemma ---"""
  shell:
    """
    cd vars/taSubInDel/stats/gemma/
    gemma -g {input.mg} -p {input.pheno} -gk 1 -o {output.relM}
    cd ../../../../
    """

rule get_vartype:
  input:
    vcf = ancient('vars/taSubInDel.{sets}.vcf')
  output:
    txt = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.typ.txt'
  message: 
    """--- Getting variant type for vcf ---"""
  shell:
    """
    grep -v '#' {input.vcf} | egrep -o 'TYP=[A-Z]+' | cut -f2 -d$'=' > {output.txt}
    """


rule sub_seg_scafbp:
  input:
    vcf = 'vars/taSubInDelBypop.vcf',
    scafbp = ancient('vars/taSubInDel/stats/gemma/vars_seg.gemma.scafbp')
  output:
    vcf = 'vars/taSubInDel.seg.vcf'
  message:
    """--- Taking subset of LD-pruned variants ---"""
  shell:
    """
    script/sub_vcf_scafbp.py {input.vcf} {input.scafbp} > {output.vcf}
    """

#rule unilmm:
#  input: 
#    geno
#    pheno
#    #envir
#    relmat
#    scafbp (-a; annot?)
#    typ
#  output:
#     directory('folder')
#  message:
#    """--- Univariate lmm assoc. test ---"""
#  shell:
#    """
#    gemma -g {input.geno} -p {input.pheno} -n 1 -gxe {input.envir} -k {input.relmat} -km 2 -lmm 4 -o {output.dir}
#    """

#rule multilmm:
#  input:
#    genoa    # -g
#    pheno    # -p 
#    envir    # -gxe
#    relmat   # -k
#    scafbp
#    #typ?
#  output:
#     directory('folder')
#  message:
#    """--- Multivariate lmm assoc. test ---"""
#  shell:
#    """
#    gemma -g  {input.geno} -p {input.pheno} -n 1 -gxe {input.envir} -lmm 4 -o {output.dir}
#    """

#rule bslmm:
#  input:
#    geno
#    pheno
#    #relmat
#    #envir
#    scafbp
#    #typ?
#  output:
#     directory('folder')
#  message:
#    """--- bslmm assoc. test ---"""
#  shell:
#    """
#    gemma -g  {input.geno} -p {input.pheno} -n 1 -gxe {input.envir} -lmm 4 -o {output.dir}
#    """

rule sub_ldp_scafbp:
  input:
    vcf = 'vars/taSubInDelBypop.vcf',
    scafbp = ancient('vars/taSubInDel/stats/al/pca/ld_prunedVars.scafbp.txt')
  output:
    vcf = 'vars/taSubInDel.ldp.vcf'
  message:
    """--- Taking subset of LD-pruned variants ---"""
  shell:
    """
    script/sub_vcf_scafbp.py {input.vcf} {input.scafbp} > {output.vcf}
    """
#

rule ngsRelate:
  input:
    vcf = 'vars/taSubInDelBypop.vcf'
  output:
    stats = 'vars/taSubInDel/stats/ngsRelate/stats.txt'
  threads: 12
  message:
    """--- Calculating relatedness (etc.) with ngsRelate2 ---"""
  shell:
    """
    module load ngsrelate/2.0
    ngsrelate -h {input.vcf} -p 12 -T GT -c 1 -O {output.stats} -z vars/ids.txt
    """

rule relStats:
  input:
    zarr = 'vars/taSubInDel.zarr/'
  output:
    touch("vars/taSubInDel/relStats.done")
  message:
    """--- Creating ngsRelate summary statistics & figures ---"""
  shell:
    """
    python script/nrel.py {input.zarr}
    """

