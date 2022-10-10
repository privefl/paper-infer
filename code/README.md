## Structure of the code

### Prepare data

- `prepare-geno-simu.R`: prepare the genotype matrix for simulations

- `prepare-geno-hm3.R`: prepare the genotype matrix for the UK Biobank analyses

- `prepare-geno-hm3-plus.R` + `greedy-maxtag.cpp`: prepare the extended SNP set proposed, HapMap3+

- `prepare-corr-hm3-*.R`: prepare the LD matrices for the three alternatives 
    - `small`: from 2000 individuals only
    - `altpop`: from individuals of around South Europe
    - `plus`: from an extended set of variants

- `prepare-pheno-fields.R` + `prepare-phecodes.R`: Prepare 248 phenotypes from fields and ICD codes in the UK Biobank


### Analyses and results

- `simu.R` + `simu-binary.R`: simulations + Figures 2 & S1--S9

- `run-GWAS.R` + `run-ldpred2-ukbb.R`: run LDpred2-auto for 248 phenotypes in the UK Biobank + Figures 3--4 & S10--S20 & S22--S23

- `run-ldpred2-cardio-proteins.R`: run LDpred2-auto for 90 external cardiovascular proteins + Figure S21

- `height-*.R`: (inference and enrichment) analyses for height with different samples sizes + Figure S24--S25

- `LPA-targeted`: predicting lipoprotein(a) concentration used a targeted region of the genome and penalized regression (individual-level data)

- `export-ldred.R`: Export precomputed LD matrices for HapMap3+
