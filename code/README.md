## Structure of the code

### Prepare data

- `prepare-geno-simu.R`: prepare the genotype matrix for simulations

- `prepare-geno-hm3.R`: prepare the genotype matrix for the UK Biobank analyses

- `prepare-geno-hm3-plus.R` + `greedy-maxtag.cpp`: prepare the extended SNP set proposed, HapMap3+

- `prepare-corr-hm3-*.R`: prepare the LD matrices for the three alternatives 
    - `small`: from 2000 individuals only
    - `altpop`: from individuals of around South Europe
    - `plus`: from an extended set of variants
    - `plus-regul`: with regularization of LD based on genetic distances

- `prepare-pheno-fields.R` + `prepare-phecodes.R`: Prepare 248 phenotypes from fields and ICD codes in the UK Biobank


### Analyses and results

- `simu-methods.R` + `simu.R` + `simu-binary.R`: simulations + Figures 2--4 & S1--S14 & S16--S23

- `compare-susie.R`: compare calibrations of LDpred2-auto and SuSiE-RSS in some simulations + Figure S15

- `run-GWAS.R` + `run-ldpred2-ukbb.R`: run LDpred2-auto for 248 phenotypes in the UK Biobank + Figures 5--6 & S24--S36 & S38--S39 & S45

- `run-ldpred2-cardio-proteins.R`: run LDpred2-auto for 90 external cardiovascular proteins + Figure S37

- `height-*.R`: (inference and enrichment) analyses for height with different samples sizes + Figures S42--S43

- `run-ldpred2-sumstats.R`: run LDpred2-auto with external GWAS summary statistics + Figure S44

- `LPA-targeted.R`: predicting lipoprotein(a) concentration used a targeted region of the genome and penalized regression (individual-level data)

- `height-one-region.R`: localized fine-mapping analysis + Figures S40--S41

- `validate-hm3-plus.R`: compute tagging statistics for HapMap3+ variants

- `export-ldred.R`: Export precomputed LD matrices for HapMap3+

- `example-with-provided-LD.R`: Example script on how to use the provided LD
