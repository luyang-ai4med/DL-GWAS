chroms = [str(i) for i in range(1,23)]
phenos = [str(i) for i in range(1,241)] # 240 phenotypes in total (80 phenotypes x prop, original, combined)
female_phenos = [str(i) for i in range(1,10)]
male_phenos = [str(i) for i in range(1,4)]

PATH_UKBB = # PATH to genotyped UKBB data
sample_list = # list of included samples
sample_list_female = # list of included female samples
sample_list_male = # list of included male samples
sample_list_5000 = # 5,000 random samples among the included samples (for predictor heritability calculation)
sample_list_20000 = # 20,000 random samples among the included samples (for LD-clumping)
snp_MAF01 = # included list of high-quality markers
annotbld65_path = # PATH to 64 BLD-LDAK annotations (http://dougspeed.com/bldldak/)

rule all:
    input:
        expand("output/Heritability/Heritability_bld65.{pheno}.hers", pheno = phenos),
        expand("output/Heritability/female_Heritability_bld65.{pheno}.hers", pheno = female_phenos),
        expand("output/Heritability/male_Heritability_bld65.{pheno}.hers", pheno = male_phenos),
        expand("output/GWAS/clump_MAF01.{pheno}.in", pheno = phenos),
        expand("output/GWAS/female_clump_MAF01.{pheno}.in", pheno = female_phenos),
        expand("output/GWAS/male_clump_MAF01.{pheno}.in", pheno = male_phenos),
        expand("output/PRS/phenotype_all_samples_PRS.{pheno}.profile", pheno = phenos),
        expand("output/PRS/female_phenotype_all_samples_PRS.{pheno}.profile", pheno = female_phenos),
        expand("output/PRS/male_phenotype_all_samples_PRS.{pheno}.profile", pheno = male_phenos),

# run a simple GWAS (linear regression) for each phenotype
rule calculate_summary_statistics:
    input:
        samples = sample_list,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        pheno = "output/phenotypes_all.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/sumstat/pheno.{pheno}"
    resources:
        mem_mb = 20000
    threads: 8
    output:       
        "output/sumstat/pheno.{pheno}.summaries"
    shell:
        "ldak5 --linear {params.out} --bfile {params.geno} --keep {input.samples} --max-threads {threads} --pheno {input.pheno} --mpheno {wildcards.pheno}"

rule clump_GWAS:
    input:
        gwas = "output/sumstat/pheno.{pheno}.pvalues",
        ref_samples = sample_list_20000,
        ref = PATH_UKBB + "ukb_snp_allChrs.bed",
        snps = snp_MAF01
    params:
        ref = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/GWAS/clump_MAF01.{pheno}"
    threads: 1
    resources:
        mem_mb = 10000
    output:
        "output/GWAS/clump_MAF01.{pheno}.in"
    shell:
        """
        ldak5 --thin-tops {params.out} --bfile {params.ref} --keep {input.ref_samples} --window-prune 0.05 --window-kb 1000  --pvalues {input.gwas} --cutoff 5e-8 --extract {input.snps}
        """

rule calculate_summary_statistics_female:
    input:
        samples = sample_list_female,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        pheno = "output/phenotypes_female.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/sumstat/female_pheno.{pheno}"
    resources:
        mem_mb = 20000
    threads: 8
    output:       
        "output/sumstat/female_pheno.{pheno}.summaries",
        "output/sumstat/female_pheno.{pheno}.pvalues"
    shell:
        "ldak5 --linear {params.out} --bfile {params.geno} --keep {input.samples} --max-threads {threads} --pheno {input.pheno} --mpheno {wildcards.pheno}"


rule clump_GWAS_female:
    input:
        gwas = "output/sumstat/female_pheno.{pheno}.pvalues",
        ref_samples = sample_list_20000,
        ref = PATH_UKBB + "ukb_snp_allChrs.bed",
        snps = snp_MAF01
    params:
        ref = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/GWAS/female_clump_MAF01.{pheno}"
    threads: 1
    resources:
        mem_mb = 10000
    output:
        "output/GWAS/female_clump_MAF01.{pheno}.in"
    shell:
        """
        ldak5 --thin-tops {params.out} --bfile {params.ref} --keep {input.ref_samples} --window-prune 0.05 --window-kb 1000  --pvalues {input.gwas} --cutoff 5e-8 --extract {input.snps}
        """

rule calculate_summary_statistics_male:
    input:
        samples = sample_list_male,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        pheno = "output/phenotypes_male.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/sumstat/male_pheno.{pheno}"
    resources:
        mem_mb = 20000
    threads: 8
    output:       
        "output/sumstat/male_pheno.{pheno}.summaries",
        "output/sumstat/male_pheno.{pheno}.pvalues",
    shell:
        "ldak5 --linear {params.out} --bfile {params.geno} --keep {input.samples} --max-threads {threads} --pheno {input.pheno} --mpheno {wildcards.pheno}"

rule clump_GWAS_male:
    input:
        gwas = "output/sumstat/male_pheno.{pheno}.pvalues",
        ref_samples = sample_list_20000,
        ref = PATH_UKBB + "ukb_snp_allChrs.bed",
        snps = snp_MAF01
    params:
        ref = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/GWAS/male_clump_MAF01.{pheno}"
    threads: 1
    resources:
        mem_mb = 10000
    output:
        "output/GWAS/male_clump_MAF01.{pheno}.in"
    shell:
        """
        ldak5 --thin-tops {params.out} --bfile {params.ref} --keep {input.ref_samples} --window-prune 0.05 --window-kb 1000  --pvalues {input.gwas} --cutoff 5e-8 --extract {input.snps}
        """

rule estimate_heritability:
    input:
        tags = "data/bld.ldak.genotyped.gbr.tagging",
        summary = "output/sumstat/pheno.{pheno}.summaries"
    params:
        out = "output/Heritability/Heritability_bld65.{pheno}"
    threads: 8
    output:   
        "output/Heritability/Heritability_bld65.{pheno}.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --check-sums NO
        """

rule estimate_heritability_female:
    input:
        tags = "data/bld.ldak.genotyped.gbr.tagging",
        summary = "output/sumstat/female_pheno.{pheno}.summaries"
    params:
        out = "output/Heritability/female_Heritability_bld65.{pheno}"
    threads: 8
    output:   
        "output/Heritability/female_Heritability_bld65.{pheno}.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --check-sums NO
        """

rule estimate_heritability_male:
    input:
        tags = "data/bld.ldak.genotyped.gbr.tagging",
        summary = "output/sumstat/male_pheno.{pheno}.summaries"
    params:
        out = "output/Heritability/male_Heritability_bld65.{pheno}"
    threads: 8
    output:   
        "output/Heritability/male_Heritability_bld65.{pheno}.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --check-sums NO
        """

# preparation to calculate per-SNP heritability (only 5000 random samples are needed)
rule estimate_prepredictor_heritabilities_preprocessing:
    input:
        samples = sample_list_5000, # 5000 are sufficient for this step
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        annot = annotbld65_path + "bld0"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        annot = annotbld65_path + "bld",
        out = "output/PRS/bld65.ldak",
        matrix = "output/PRS/bld65.ldak.matrix",
    threads: 8
    output:       
        weights = annotbld65_path + "bld65"
    shell:
        """
        ldak5 --cut-weights sections --bfile {params.geno} --keep {input.samples} --max-threads {threads}
        ldak5 --calc-weights-all sections --bfile {params.geno} --keep {input.samples} --max-threads {threads}
        mv sections/weights.short {output.weights}        
        """
# calculate tagging files based on 65 annotations (same 5000 samples)
rule calculate_snp_tags:
    input:
        samples = sample_list_5000, # 5000 are sufficient for this step
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        annot = annotbld65_path + "bld0",
        weights = annotbld65_path + "bld65",
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        annot = annotbld65_path + "bld",
        out = "output/PRS/bld65.ldak"
    threads: 8
    output:   
        tags = "output/PRS/bld65.ldak.tagging",
        matrix = "output/PRS/bld65.ldak.matrix"
    shell:
        """
        ldak5 --calc-tagging {params.out} --bfile {params.geno} --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix {params.annot} --window-kb 1200 --save-matrix YES --keep {input.samples} --max-threads {threads}
        """

# calculate per-predictor heritabilities which will be the input for PRS estimation
# This step also calculates heritabilities for the 5000 samples.
rule estimate_prepredictor_heritabilities:
    input:
        tags = "output/PRS/bld65.ldak.tagging",
        hermat = "output/PRS/bld65.ldak.matrix",
        summary = "output/sumstat/pheno.{pheno}.summaries"
    params:
        out = "output/PRS/bld65.ldak.{pheno}"
    threads: 8
    output:   
        pred_her = "output/PRS/bld65.ldak.{pheno}.ind.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --matrix {input.hermat}
        """

# Calculate PRS effect estimates (i.e. training of the PRS model using cross-validation)
rule estimate_PRS_predictors:
    input:
        pheno = "output/phenotypes_all.txt",
        samples = sample_list,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        indhers = "output/PRS/bld65.ldak.{pheno}.ind.hers"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 10000
    output:
        "output/PRS/phenotype_all_samples_PRS.{pheno}.effects"
    shell:
        """
        ldak5 --bolt {params.out} --pheno {input.pheno} --mpheno {wildcards.pheno} --bfile {params.geno} --ind-hers {input.indhers} --cv-proportion 0.1 --keep {input.samples} --max-threads {threads}
        """

rule estimate_prepredictor_heritabilities_female:
    input:
        tags = "output/PRS/bld65.ldak.tagging",
        hermat = "output/PRS/bld65.ldak.matrix",
        summary = "output/sumstat/female_pheno.{pheno}.summaries"
    params:
        out = "output/PRS/female_bld65.ldak.{pheno}"
    threads: 8
    output:   
        pred_her = "output/PRS/female_bld65.ldak.{pheno}.ind.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --matrix {input.hermat} --check-sums NO
        """

rule estimate_prepredictor_heritabilities_male:
    input:
        tags = "output/PRS/bld65.ldak.tagging",
        hermat = "output/PRS/bld65.ldak.matrix",
        summary = "output/sumstat/male_pheno.{pheno}.summaries"
    params:
        out = "output/PRS/male_bld65.ldak.{pheno}"
    threads: 8
    output:   
        pred_her = "output/PRS/male_bld65.ldak.{pheno}.ind.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --matrix {input.hermat} --check-sums NO
        """

rule estimate_PRS_predictors_female:
    input:
        pheno = "output/phenotypes_female.txt",
        samples = sample_list_female,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        indhers = "output/PRS/female_bld65.ldak.{pheno}.ind.hers"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/female_phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 25000
    output:
        "output/PRS/female_phenotype_all_samples_PRS.{pheno}.effects"
    shell:
        """
        ldak5 --bolt {params.out} --pheno {input.pheno} --mpheno {wildcards.pheno} --bfile {params.geno} --ind-hers {input.indhers} --cv-proportion 0.1 --keep {input.samples} --max-threads {threads}
        """

rule estimate_PRS_predictors_male:
    input:
        pheno = "output/phenotypes_male.txt",
        samples = sample_list_male,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        indhers = "output/PRS/male_bld65.ldak.{pheno}.ind.hers"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/male_phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 20000
    output:
        "output/PRS/male_phenotype_all_samples_PRS.{pheno}.effects"
    shell:
        """
        ldak5 --bolt {params.out} --pheno {input.pheno} --mpheno {wildcards.pheno} --bfile {params.geno} --ind-hers {input.indhers} --cv-proportion 0.1 --keep {input.samples} --max-threads {threads}
        """

rule calculate_validate_PRS_scores:
    input:
        estimators = "output/PRS/phenotype_all_samples_PRS.{pheno}.effects",
        samples = sample_list,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        test = "output/phenotypes_future_test.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 10000
    output:
        "output/PRS/phenotype_all_samples_PRS.{pheno}.profile"
    shell:
        """
        ldak5 --calc-scores {params.out} --scorefile {input.estimators} --bfile {params.geno} --keep {input.samples} --power 0 --max-threads {threads} --pheno {input.test} --mpheno {wildcards.pheno}
        """
    
rule calculate_validate_PRS_scores_female:
    input:
        estimators = "output/PRS/female_phenotype_all_samples_PRS.{pheno}.effects",
        samples = sample_list_female,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        test = "output/phenotypes_future_test.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/female_phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 10000
    output:
        "output/PRS/female_phenotype_all_samples_PRS.{pheno}.profile"
    shell:
        """
        ldak5 --calc-scores {params.out} --scorefile {input.estimators} --bfile {params.geno} --keep {input.samples} --power 0 --max-threads {threads} --pheno {input.test} --mpheno {wildcards.pheno}
        """

rule calculate_validate_PRS_scores_male:
    input:
        estimators = "output/PRS/male_phenotype_all_samples_PRS.{pheno}.effects",
        samples = sample_list_male,
        geno = PATH_UKBB + "ukb_snp_allChrs.bed",
        test = "output/phenotypes_future_test.txt"
    params:
        geno = PATH_UKBB + "ukb_snp_allChrs",
        out = "output/PRS/male_phenotype_all_samples_PRS.{pheno}"
    threads: 8
    resources:
        mem_mb = 10000
    output:
        "output/PRS/male_phenotype_all_samples_PRS.{pheno}.profile",
    shell:
        """
        ldak5 --calc-scores {params.out} --scorefile {input.estimators} --bfile {params.geno} --keep {input.samples} --power 0 --max-threads {threads} --pheno {input.test} --mpheno {wildcards.pheno}
        """