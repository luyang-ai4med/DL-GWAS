import pandas as pd

chroms = [str(i) for i in range(1,23)]

PATH_PHENO = "/oak/stanford/groups/rbaltman/luyang/scripts/deepGWAS/revision/results_full3/"
pheno_map_file = "/oak/stanford/groups/rbaltman/luyang/scripts/deepGWAS/revision/popdx_train_phenotypes.csv"
cov_file = "/oak/stanford/groups/rbaltman/ukbiobank/phenotype_data/gwas_files/ukb33722_GWAS_covar.phe"
pheno_df = pd.read_csv(pheno_map_file)
phenos = list(pheno_df["index"].astype(str))

PATH_UKBB = "/oak/stanford/groups/rbaltman/ukbiobank/genotype_data/all_genotype_files/"
sample_list = "/oak/stanford/groups/rbaltman/luyang/scripts/Heritability_v2/data/UKBB_British_unrelated_samples_phenotypes_all.txt"
sample_list_female = "data/UKBB_British_unrelated_samples_phenotypes_female.txt"
sample_list_male = "data/UKBB_British_unrelated_samples_phenotypes_male.txt"
sample_list_5000 = "data/UKBB_British_unrelated_samples_phenotypes_all_5000.txt"
sample_list_20000 = "data/random_samples_20000.txt"

snp_QC = "output_all_features/snp_qc.snplist"
annotbld65_path = "/oak/stanford/groups/rbaltman/luyang/scripts/Heritability_v2/data/bld/"
tagging_file = "data/bld.ldak.genotyped.gbr.tagging"
prevalence_file = "data/popdx_train_phenotypes_prevalence.csv"

rule all:
    input:
        expand("output_all_features/pheno/"),
        expand("output_all_features/Heritability/Heritability_bld65.{pheno}.{set}.{type}.hers", pheno = phenos, set = ["all", "train"], type = ["quantitative", "binary", "meta"]),
        expand("output_all_features/Heritability/Heritability_bld65.{pheno}.{set}.{type}.hers", pheno = phenos, set = ["all"], type = ["quantitativeall"]),
        expand("output_all_features/GWAS/pheno.{pheno}.{set}.{type}.in", pheno = phenos, set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        #expand("output_all_features/LD/pheno.{pheno}.{set}.{type}.ld", pheno = phenos, set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/LD/pheno.{pheno}.{set}.{type}.0p6.ld", pheno = phenos, set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/PRS/sample_all.{pheno}.{set}.{type}.profile", pheno = phenos, set = ["all"], type =["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/PRS/sample_all.{pheno}.{set}.{type}.profile", pheno = phenos, set = ["train"], type = ["quantitative", "binary", "meta"]),
        expand("output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitativeall_binary.cors", pheno = phenos, set = ["all"]),
        expand("output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitative_binary.cors", pheno = phenos, set = ["all"]),
        expand("output_all_features/ldsc/LDSC.{pheno}.{set}.quantitativeall_binary.log", pheno = phenos, set = ["all"]),
        expand("output_all_features/ldsc/LDSC.{pheno}.{set}.quantitative_binary.log", pheno = phenos, set = ["all"]),
        expand("output_all_features/ldsc/LDSC.{pheno}.{set}.meta_binary.log", pheno = phenos, set = ["all"]),
        expand("output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.summaries", pheno = phenos, set = ["all"]),
        expand("output_all_features/ldsc/spa.LDSC.{pheno}.{set}.quantitativeall.spa.log", pheno = phenos, set = ["all"]),
        #expand("output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitativeall_edu.log", pheno = phenos, set = ["all"]),
        #expand("output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitative_edu.log", pheno = phenos, set = ["all"]),
        #expand("output_all_features/edu_ldsc/LDSC.{pheno}.{set}.binary_edu.log", pheno = phenos, set = ["all"]),
        #expand("output_all_features/edu_ldsc/LDSC.{pheno}.{set}.meta_edu.log", pheno = phenos, set = ["all"]),
        expand("output_all_features/Heritability/Heritability_binary_ascertainment.csv"),
        expand("output_all_features/validation_ldsc/LDSC.672.{set}.{type}_AF.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.780.{set}.{type}_asthma.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.308.{set}.{type}_obesity.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.92.{set}.{type}_bc.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.632.{set}.{type}_cad.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.532.{set}.{type}_cataract.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.275.{set}.{type}_ldl.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.1316.{set}.{type}_oa.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.623.{set}.{type}_sbp.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.697.{set}.{type}_stroke.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"]),
        expand("output_all_features/validation_ldsc/LDSC.205.{set}.{type}_t2d.log", set = ["all"], type = ["quantitativeall","quantitative", "binary", "meta"])
        

rule prepare_phenotype:
    input:
        cov = cov_file,
        samples = sample_list,
        pheno_file = pheno_map_file,
    params:
        PATH = PATH_PHENO
    resources:
        mem_mb = 5000
    output:
        all_q = "output_all_features/pheno/{pheno}.all.quantitative",
        all_all_q = "output_all_features/pheno/{pheno}.all.quantitativeall",
        all_b = "output_all_features/pheno/{pheno}.all.binary",
        all_cov = "output_all_features/pheno/{pheno}.all.cov",
        train_q = "output_all_features/pheno/{pheno}.train.quantitative",
        train_b = "output_all_features/pheno/{pheno}.train.binary",
        train_cov = "output_all_features/pheno/{pheno}.train.cov",
    script:
        "scripts/prepare_phenotype_all.R"

# run a simple GWAS (linear regression) for controls with disease liabilities
rule calculate_quantitative_summary_statistics:
    input:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed", 
        pheno = "output_all_features/pheno/{pheno}.{set}.quantitative"
    params:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitative"
    resources:
        mem_mb = 15000
    threads: 8
    output:       
        "output_all_features/sumstat/pheno.{pheno}.{set}.quantitative.summaries"
    shell:
        "ldak5 --linear {params.out} --bfile {params.geno} --max-threads {threads} --pheno {input.pheno}"

# run a simple GWAS (linear regression) for controls & cases with disease liabilities
rule calculate_quantitative_summary_statistics_all:
    input:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed", 
        pheno = "output_all_features/pheno/{pheno}.{set}.quantitativeall"
    params:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall"
    resources:
        mem_mb = 15000
    threads: 8
    output:       
        "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.summaries"
    shell:
        "ldak5 --linear {params.out} --bfile {params.geno} --max-threads {threads} --pheno {input.pheno}"




# run a binary GWAS 
rule calculate_binary_summary_statistics_ldak:
    input:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed", 
        pheno = "output_all_features/pheno/{pheno}.{set}.binary",
        cov = "output_all_features/pheno/{pheno}.{set}.cov"
    params:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/sumstat/pheno.{pheno}.{set}.binary"
    resources:
        mem_mb = 15000
    threads: 8
    output:       
        "output_all_features/sumstat/pheno.{pheno}.{set}.binary.summaries"
    shell:
        "ldak5 --logistic {params.out} --bfile {params.geno} --max-threads {threads} --pheno {input.pheno} --covar {input.cov}"



### meta-analysis of quantitative GWAS on control and binary GWAS
rule meta_analyse_GWAS:
    input:
        quant = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitative.summaries",
        binary = "output_all_features/sumstat/pheno.{pheno}.{set}.binary.summaries"
    resources:
        mem_mb = 10000
    threads: 1
    output:
        "output_all_features/sumstat/pheno.{pheno}.{set}.meta.summaries",
        "output_all_features/sumstat/pheno.{pheno}.{set}.meta.pvalues"
    script:
        "scripts/meta_analyse_GWAS.R"


## clump GWAS results based on thresholds
rule clump_GWAS:
    input:
        pvals = "output_all_features/sumstat/pheno.{pheno}.{set}.{type}.pvalues",
        ref_samples = sample_list_20000,
        ref = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed",
        snps = snp_QC
    params:
        ref = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/GWAS/pheno.{pheno}.{set}.{type}"
    threads: 1
    resources:
        mem_mb = 10000
    output:
        temp_out = temp("output_all_features/GWAS/pheno.{pheno}.{set}.{type}.tempsig"),
        temp_out2 = temp("output_all_features/GWAS/pheno.{pheno}.{set}.{type}.tempsigfilter"),
        out = "output_all_features/GWAS/pheno.{pheno}.{set}.{type}.in"
    shell:
        """
        awk '{{if ($2 <= 5e-8) {{ print $1 }} }}' {input.pvals} > {output.temp_out}
        awk 'FNR==NR{{col1[$1]; next}}($1 in col1){{print $1}}' {output.temp_out} {input.snps} > {output.temp_out2}
        if [ -s {output.temp_out2} ]
        then
            ldak5 --thin-tops {params.out} --bfile {params.ref} --keep {input.ref_samples} --window-prune 0.05 --window-kb 1000  --pvalues {input.pvals} --cutoff 5e-8 --extract {input.snps}
        else
            touch {output.out}
        fi  
        """


## LD snps  
rule LD_clump:
    input:
        clump_snps = "output_all_features/GWAS/pheno.{pheno}.{set}.{type}.in",
        ref_samples = sample_list_20000,
        ref = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed",
        snps = snp_QC
    params:
        ref = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/LD/pheno.{pheno}.{set}.{type}.0p6"
    threads: 1
    resources:
        mem_mb = 10000
    output:
        out = "output_all_features/LD/pheno.{pheno}.{set}.{type}.0p6.ld"
    shell:
        """
        if [ -s {input.clump_snps} ]
        then
            plink --bfile {params.ref} --keep {input.ref_samples} --r2 --ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 0.6 \
            --out {params.out} --ld-snp-list {input.clump_snps}
        else
            touch {output.out}
        fi
        """

rule estimate_heritability:
    input:
        tags = tagging_file,
        summary = "output_all_features/sumstat/pheno.{pheno}.{set}.{type}.summaries",
        snps = snp_QC
    params:
        out = "output_all_features/Heritability/Heritability_bld65.{pheno}.{set}.{type}"
    threads: 8
    output:   
        "output_all_features/Heritability/Heritability_bld65.{pheno}.{set}.{type}.hers"
    shell:
        """
        ldak5 --sum-hers {params.out} --tagfile {input.tags} --summary {input.summary} --extract {input.snps} --check-sums NO 
        """

#prepare ascertainment (proportion of cases) to convert observed scale h2 to liability scale h2
rule calculate_binary_ascertainment:
    input:
        pheno_file = pheno_map_file,
        pheno_path = "output_all_features/pheno/",
        prevalence_file = prevalence_file
    output:
        "output_all_features/Heritability/Heritability_binary_ascertainment.csv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/calculate_ascertainment.py"

rule genetic_correlation_binary_quantitativeall_ldak:
    input:
        tags = tagging_file,
        summary1 = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.summaries",
        summary2 = "output_all_features/sumstat/pheno.{pheno}.{set}.binary.summaries",
        snps = snp_QC
    params:
        out = "output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitativeall_binary"
    threads: 1
    resources:
        mem_mb = 5000
    output:   
        "output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitativeall_binary.cors"
    shell:
        """
        ldak5 --sum-cors {params.out} --tagfile {input.tags} --summary {input.summary1} --summary2 {input.summary2} --check-sums NO --cutoff 0.01 --extract {input.snps}
        """
        
        
        
rule genetic_correlation_binary_quantitative_ldak:
    input:
        tags = tagging_file,
        summary1 = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitative.summaries",
        summary2 = "output_all_features/sumstat/pheno.{pheno}.{set}.binary.summaries",
        snps = snp_QC
    params:
        out = "output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitative_binary"
    threads: 1
    resources:
        mem_mb = 5000
    output:   
        "output_all_features/Correlation/Correlation_bld65.{pheno}.{set}.quantitative_binary.cors"
    shell:
        """
        ldak5 --sum-cors {params.out} --tagfile {input.tags} --summary {input.summary1} --summary2 {input.summary2} --check-sums NO --cutoff 0.01 --extract {input.snps}
        """


# calculate per-predictor heritabilities which will be the input for PRS estimation with pre-computed tagging file
rule calculate_prepredictor_heritabilities:
    input:
        sumstat = "output_all_features/sumstat/pheno.{pheno}.{set}.{type}.summaries",    ##train 
        tagging = "data/gbr.genotyped/gbr.genotyped.bld.ldak.quickprs.tagging",
        matrix = "data/gbr.genotyped/gbr.genotyped.bld.ldak.quickprs.matrix",
        snps = snp_QC,

    params:
        "output_all_features/PRS/pheno.{pheno}.{set}.{type}"
    resources:
        mem_mb = 5000
    threads: 4
    output:
        "output_all_features/PRS/pheno.{pheno}.{set}.{type}.ind.hers"
    shell:
        """
        ldak5 --sum-hers {params} --summary {input.sumstat} --tagfile {input.tagging} \
        --matrix {input.matrix} --check-sums NO --max-threads {threads} --extract {input.snps}
        """

# estimate PRS effects (QUICK-PRS method)
rule calculate_PRS_effects:
    input:
        sumstat = "output_all_features/sumstat/pheno.{pheno}.{set}.{type}.summaries", ##train
        her = "output_all_features/PRS/pheno.{pheno}.{set}.{type}.ind.hers",
        cor = "data/gbr.genotyped/gbr.genotyped.cors.bin",
        ld = "data/gbr.genotyped/highld.snps",
        snpexclude = "data/snp.exclude",
        snps = snp_QC
    params:
        cor = "data/gbr.genotyped/gbr.genotyped",
        out = "output_all_features/PRS/pheno.{pheno}.{set}.{type}.bayesr"
    resources:
        mem_mb = 10000
    threads: 4
    output:
        "output_all_features/PRS/pheno.{pheno}.{set}.{type}.bayesr.effects"
    shell:
        """
        ldak5 --mega-prs {params.out} --summary {input.sumstat} --ind-hers {input.her} --cors {params.cor} \
        --high-LD {input.ld} --model bayesr --cv-proportion .1 --window-cm 1 --max-threads {threads} \
        --exclude {input.snpexclude} --extract {input.snps}
        """

rule calculate_PRS_scores:
    input:
        estimators = "output_all_features/PRS/pheno.{pheno}.{set}.{type}.bayesr.effects",
        samples = sample_list,
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed",
        snps = snp_QC,
    params:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/PRS/sample_all.{pheno}.{set}.{type}"
    threads: 8
    resources:
        mem_mb = 10000
    output:
        "output_all_features/PRS/sample_all.{pheno}.{set}.{type}.profile"
    shell:
        """
        ldak5 --calc-scores {params.out} --scorefile {input.estimators} --bfile {params.geno} --keep {input.samples} --power 0 --max-threads {threads} --extract {input.snps}
        """

## run ldsc correlations
rule format_munge:
    input:
        GWAS = "output_all_features/sumstat/pheno.{pheno}.{set}.{type}.summaries",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/ldsc/pheno.{pheno}.{set}.{type}.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge.R"

rule munge_sumstat:
    input:
        GWAS = "output_all_features/ldsc/pheno.{pheno}.{set}.{type}.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/ldsc/pheno.{pheno}.{set}.{type}.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/ldsc/pheno.{pheno}.{set}.{type}.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_genetic_correlation_quantitativeall_binary:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.quantitativeall.munge.sumstats.gz",
        GWAS2 = "output_all_features/ldsc/pheno.{pheno}.{set}.binary.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/ldsc/LDSC.{pheno}.{set}.quantitativeall_binary"
    conda: "ldsc"
    output:
        "output_all_features/ldsc/LDSC.{pheno}.{set}.quantitativeall_binary.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"
        
rule ldsc_genetic_correlation_quantitative_binary:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.quantitative.munge.sumstats.gz",
        GWAS2 = "output_all_features/ldsc/pheno.{pheno}.{set}.binary.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/ldsc/LDSC.{pheno}.{set}.quantitative_binary"
    conda: "ldsc"
    output:
        "output_all_features/ldsc/LDSC.{pheno}.{set}.quantitative_binary.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"
        
rule ldsc_genetic_correlation_meta_binary:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.meta.munge.sumstats.gz",
        GWAS2 = "output_all_features/ldsc/pheno.{pheno}.{set}.binary.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/ldsc/LDSC.{pheno}.{set}.meta_binary"
    conda: "ldsc"
    output:
        "output_all_features/ldsc/LDSC.{pheno}.{set}.meta_binary.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"  
        
        
#ldak spa     
rule calculate_quantitative_summary_statistics_all_spa:
    input:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs.bed", 
        pheno = "output_all_features/pheno/{pheno}.{set}.quantitativeall"
    params:
        geno = f'{PATH_UKBB}' + "ukb_snp_allChrs",
        out = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear"
    resources:
        mem_mb = 15000
    threads: 8
    output:       
        "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.summaries",
        "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.spa"
    shell:
        "ldak5.2.linux --linear {params.out} --bfile {params.geno} --max-threads {threads} --pheno {input.pheno} --spa-test YES"
        
## run ldsc correlations for spa and non-spa quantitativeall
rule format_munge_quantitativeall:
    input:
        GWAS = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.summaries",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge.R"
        
rule qc_and_format_spa_GWAS:
    input:
        gwas = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.spa",
    resources:
        mem_mb = 10000
    threads: 1
    output:
        "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.spa.summaries",
    script:
        "scripts/qc_and_format_spa.py"
        
rule format_munge_quantitativeall_spa:
    input:
        GWAS = "output_all_features/sumstat/pheno.{pheno}.{set}.quantitativeall.linear.spa.summaries",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.spa.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge.R"

rule munge_sumstat_quantitativeall:
    input:
        GWAS = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"
        
rule munge_sumstat_quantitativeall_spa:
    input:
        GWAS = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.spa.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.spa.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.spa.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_genetic_correlation_quantitativeall_spa:
    input:
        GWAS1 = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.munge.sumstats.gz",
        GWAS2 = "output_all_features/ldsc/spa.pheno.{pheno}.{set}.quantitativeall.linear.spa.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/ldsc/spa.LDSC.{pheno}.{set}.quantitativeall.spa"
    conda: "ldsc"
    output:
        "output_all_features/ldsc/spa.LDSC.{pheno}.{set}.quantitativeall.spa.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


## ldak education 
        
rule format_education:
    input:
        GWAS_edu = "data/22501_irnt.gwas.imputed_v3.both_sexes_variants.tsv",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/edu_ldsc/pheno.edu.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_education.R"
        
rule education_sumstat:
    input:
        GWAS = "output_all_features/edu_ldsc/pheno.edu.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/edu_ldsc/pheno.edu.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/edu_ldsc/pheno.edu.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_quantitativeall_edu:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.quantitativeall.munge.sumstats.gz",
        GWAS2 = "output_all_features/edu_ldsc/pheno.edu.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitativeall_edu"
    conda: "ldsc"
    output:
        "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitativeall_edu.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"
        
rule ldsc_quantitative_edu:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.quantitative.munge.sumstats.gz",
        GWAS2 = "output_all_features/edu_ldsc/pheno.edu.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitative_edu"
    conda: "ldsc"
    output:
        "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.quantitative_edu.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"
        
rule ldsc_binary_edu:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.binary.munge.sumstats.gz",
        GWAS2 = "output_all_features/edu_ldsc/pheno.edu.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.binary_edu"
    conda: "ldsc"
    output:
        "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.binary_edu.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"      
        
rule ldsc_meta_edu:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.{pheno}.{set}.meta.munge.sumstats.gz",
        GWAS2 = "output_all_features/edu_ldsc/pheno.edu.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.meta_edu"
    conda: "ldsc"
    output:
        "output_all_features/edu_ldsc/LDSC.{pheno}.{set}.meta_edu.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"     


####### validation
#AF

rule format_AF:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/AF_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.AF.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule AF_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.AF.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.AF.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.AF.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_AF:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.672.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.AF.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.672.{set}.{type}_AF"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.672.{set}.{type}_AF.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Asthma 

rule format_asthma:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/ASTHMA_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.asthma.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule asthma_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.asthma.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.asthma.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.asthma.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_asthma:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.780.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.asthma.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.780.{set}.{type}_asthma"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.780.{set}.{type}_asthma.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"

#Obesity

rule format_obesity:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/BMI_meta_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.obesity.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule obesity_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.obesity.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.obesity.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.obesity.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_obesity:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.308.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.obesity.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.308.{set}.{type}_obesity"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.308.{set}.{type}_obesity.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"

#Breast cancer

rule format_bc:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/BREASTCANCER_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.bc.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule bc_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.bc.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.bc.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.bc.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_bc:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.92.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.bc.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.92.{set}.{type}_bc"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.92.{set}.{type}_bc.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Coronary atherosclerosis

rule format_cad:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/CAD_consortium_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.cad.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule cad_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.cad.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.cad.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.cad.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_cad:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.632.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.cad.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.632.{set}.{type}_cad"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.632.{set}.{type}_cad.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#CATARACT

rule format_cataract:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/CATARACT_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.cataract.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule cataract_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.cataract.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.cataract.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.cataract.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_cataract:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.532.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.cataract.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.532.{set}.{type}_cataract"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.532.{set}.{type}_cataract.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Hypercholesterolemia

rule format_ldl:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/LDL_EUR_2021_summary_stats.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.ldl.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule ldl_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.ldl.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.ldl.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.ldl.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_ldl:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.275.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.ldl.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.275.{set}.{type}_ldl"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.275.{set}.{type}_ldl.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Osteoarthritis

rule format_oa:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/OA_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.oa.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule oa_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.oa.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.oa.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.oa.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_oa:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.1316.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.oa.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.1316.{set}.{type}_oa"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.1316.{set}.{type}_oa.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Essential hypertension

rule format_sbp:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/SBP_meta_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.sbp.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule sbp_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.sbp.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.sbp.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.sbp.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_sbp:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.623.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.sbp.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.623.{set}.{type}_sbp"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.623.{set}.{type}_sbp.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"


#Stroke

rule format_stroke:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/STROKE_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.stroke.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule stroke_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.stroke.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.stroke.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.stroke.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_stroke:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.697.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.stroke.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.697.{set}.{type}_stroke"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.697.{set}.{type}_stroke.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"



#Type 2 diabetes

rule format_t2d:
    input:
        GWAS_external = "/oak/stanford/groups/rbaltman/luyang/scripts/popdx_genetics_revision/Heritability_v2_all_features/data/GWAS/T2D_gwas_summary_uk10kck.ma",
        snp_map = f'{PATH_UKBB}' + "ukb_snp_allChrs_snp_mapping.tsv"
    output:
        "output_all_features/validation_ldsc/pheno.t2d.munge.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/format_munge_external.R"
        
rule t2d_sumstat:
    input:
        GWAS = "output_all_features/validation_ldsc/pheno.t2d.munge.tsv",
        snps = "data/ldsc/w_hm3.snplist"
    params:
        out = "output_all_features/validation_ldsc/pheno.t2d.munge"
    threads: 1
    resources:
        mem_mb = 5000
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/pheno.t2d.munge.sumstats.gz"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/munge_sumstats.py --chunksize 500000 --sumstats {input.GWAS} --out {params.out} --merge-alleles {input.snps}"

rule ldsc_t2d:
    input:
        GWAS1 = "output_all_features/ldsc/pheno.205.{set}.{type}.munge.sumstats.gz",
        GWAS2 = "output_all_features/validation_ldsc/pheno.t2d.munge.sumstats.gz",
        ld_scores = "data/ldsc/eur_w_ld_chr/1.l2.M_5_50"
    params:
        ld_scores = "data/ldsc/eur_w_ld_chr/",
        out = "output_all_features/validation_ldsc/LDSC.205.{set}.{type}_t2d"
    conda: "ldsc"
    output:
        "output_all_features/validation_ldsc/LDSC.205.{set}.{type}_t2d.log"
    shell:
        "/oak/stanford/groups/rbaltman/luyang/software/ldsc/ldsc.py --rg {input.GWAS1},{input.GWAS2} --ref-ld-chr {params.ld_scores} --w-ld-chr {params.ld_scores} --out {params.out}"

