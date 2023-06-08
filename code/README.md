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

- `simu-methods.R` + `simu.R` + `simu-binary.R`: simulations + Figures 2 & S1--S7 & S9--S14

- `compare-susie.R`: compare calibrations of LDpred2-auto and SuSiE-RSS in some simulations + Figure S8

- `run-GWAS.R` + `run-ldpred2-ukbb.R`: run LDpred2-auto for 248 phenotypes in the UK Biobank + Figures 3--4 & S15--S26 & S28--S29

- `run-ldpred2-cardio-proteins.R`: run LDpred2-auto for 90 external cardiovascular proteins + Figure S27

- `height-*.R`: (inference and enrichment) analyses for height with different samples sizes + Figure S30--S31

- `run-ldpred2-sumstats.R`: run LDpred2-auto with external GWAS summary statistics + Figure S32

- `LPA-targeted`: predicting lipoprotein(a) concentration used a targeted region of the genome and penalized regression (individual-level data)

- `export-ldred.R`: Export precomputed LD matrices for HapMap3+

- `example-with-provided-LD.R`: Example script on how to use the provided LD
