#### From https://doi.org/10.1016/j.ajhg.2021.11.008 ####

library(bigreadr)

csv <- "UKBB/ukb41181.csv"

sex <- fread2(csv, select = "22001-0.0")[[1]]


#### Quantitative UKBB phenotypes ####

bigassertr::assert_dir("data/ukbb-quant-pheno")

#### Misc ####

# Hand grip strength (mean of both) (46/47)
y <- rowMeans(fread2(csv, select = c("46-0.0", "47-0.0")))
y[y < 1] <- NA
saveRDS(y, "data/ukbb-quant-pheno/hand_grip_strength.rds")

# Waist circumference (48)
y <- fread2(csv, select = "48-0.0")[[1]]
y[log(y) < 3.5] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_waist_circ.rds")

# Hip circumference (49)
y <- fread2(csv, select = "49-0.0")[[1]]
y[log(y) < 4] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_hip_circ.rds")

# Standing height (50)
y <- fread2(csv, select = "50-0.0")[[1]]
y[y < 130] <- NA
saveRDS(y, "data/ukbb-quant-pheno/height.rds")

# sex specific
yM <- y; yM[sex == 0] <- NA
saveRDS(yM, "data/ukbb-quant-pheno/M_height.rds")
yF <- y; yF[sex == 1] <- NA
saveRDS(yF, "data/ukbb-quant-pheno/F_height.rds")

# Pulse rate, automated reading (102)
y <- fread2(csv, select = "102-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_pulse_rate.rds")

# Average total household income before tax (738)
y <- fread2(csv, select = "738-0.0")[[1]]
y[y < 0] <- NA
saveRDS(sqrt(y), "data/ukbb-quant-pheno/income.rds")

# Time spent watching television (TV) or using computer (1070/1080)
y <- fread2(csv, select = c("1070-0.0", "1080-0.0"))
y[y == -10] <- 0; y[y < 0] <- NA
y <- rowSums(y)
y[y > 15] <- NA
saveRDS(y, "data/ukbb-quant-pheno/geek_time.rds")

# Sleep duration (1160)
y <- fread2(csv, select = "1160-0.0")[[1]]
y[y < 4 | y > 12] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_sleep.rds")

# Morning/evening person (chronotype) (1180)
y <- fread2(csv, select = "1180-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/more_evening.rds")

# Sleeplessness / insomnia (1200)
y <- fread2(csv, select = "1200-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/insomnia.rds")

# Water intake (1528)
y <- fread2(csv, select = "1528-0.0")[[1]]
y[y == -10] <- 0; y[y < 0 | y > 15] <- NA
saveRDS(y, "data/ukbb-quant-pheno/water_intake.rds")

# Alcohol intake frequency (1558)
y <- fread2(csv, select = "1558-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/less_alcohol.rds")

# Skin colour (1717)
y <- fread2(csv, select = "1717-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/darker_skin.rds")

# Ease of skin tanning (1727)
y <- fread2(csv, select = "1727-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/less_tanned.rds")

# Hair colour (natural, before greying) (1747)
y <- fread2(csv, select = "1747-0.0")[[1]]
y <- dplyr::case_when(y == 1 ~ 1L,
                      y == 3 ~ 2L,
                      y == 4 ~ 3L,
                      y == 5 ~ 4L,
                      TRUE ~ NA_integer_)
saveRDS(y, "data/ukbb-quant-pheno/darker_hair.rds")

# Age first had sexual intercourse (2139)
y <- fread2(csv, select = "2139-0.0")[[1]]
y[y < 10 | y > 40] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_age_first_sex.rds")

# Overall health rating (2178)
y <- fread2(csv, select = "2178-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/poorer_health.rds")

# Use of sun/uv protection (2267)
y <- fread2(csv, select = "2267-0.0")[[1]]
y[y < 0 | y == 5] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/more_sunscreen.rds")

# Hair/balding pattern (2395)
y <- fread2(csv, select = "2395-0.0")[[1]]
y[y < 0 | sex == 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/M_less_hair.rds")

# Age when periods started (menarche) (2714)
y <- fread2(csv, select = "2714-0.0")[[1]]
y[y < 0 | sex == 1] <- NA
saveRDS(y, "data/ukbb-quant-pheno/F_menarche.rds")

# Age at first live birth (2754)
y <- fread2(csv, select = "2754-0.0")[[1]]
y[y < 0 | sex == 1] <- NA
saveRDS(y, "data/ukbb-quant-pheno/F_age_first_birth.rds")

# Ankle spacing width (3143)
y <- fread2(csv, select = "3143-0.0")[[1]]
y[log(y) < 3.2] <- NA
saveRDS(y, "data/ukbb-quant-pheno/ankle_spacing.rds")

# Heel Broadband ultrasound attenuation, direct entry (3144)
y <- fread2(csv, select = "3144-0.0")[[1]]
y[log(y) < 3] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_heel_BUA.rds")

# Speed of sound through heel (3146)
y <- fread2(csv, select = "3146-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_heel_SoS.rds")

# Heel bone mineral density (BMD) (3148)
y <- fread2(csv, select = "3148-0.0")[[1]]
y[log(y) < -2] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_heel_BMD.rds")

# Length of menstrual cycle (3710)
y <- fread2(csv, select = "3710-0.0")[[1]]
y[y < 0 | sex == 1] <- NA
y[y < 15 | y > 40] <- NA
saveRDS(y, "data/ukbb-quant-pheno/F_length_menstrual_cycle.rds")

# Diastolic blood pressure, automated reading (4079)
y <- fread2(csv, select = "4079-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/diastolic_BP.rds")

# Systolic blood pressure, automated reading (4080)
y <- fread2(csv, select = "4080-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/systolic_BP.rds")

# logMAR in round (left/right) (5078)
y <- fread2(csv, select = c(paste0("5078-0.", 0:15),
                            paste0("5078-1.", 0:15),
                            paste0("5079-0.", 0:15),
                            paste0("5079-1.", 0:15)))
y <- rowMeans(y, na.rm = TRUE)
saveRDS(y, "data/ukbb-quant-pheno/logMAR.rds")

# Qualifications (6138)
y <- fread2(csv, select = "6138-0.0")[[1]]
y <- dplyr::case_when(y == 3 ~ 10,
                      y == 4 ~ 10,
                      y == 2 ~ 13,
                      y == 6 ~ 15,
                      y == 5 ~ 19,
                      y == 1 ~ 20)
saveRDS(y, "data/ukbb-quant-pheno/years_of_edu.rds")

# Sitting height (20015)
y <- fread2(csv, select = "20015-0.0")[[1]]
y[y < 70 | y > 110] <- NA
saveRDS(y, "data/ukbb-quant-pheno/sitting_height.rds")

# Birth weight (20022)
y <- fread2(csv, select = "20022-0.0")[[1]]
y[y < 1 | y > 6] <- NA
saveRDS(y, "data/ukbb-quant-pheno/birth_weight.rds")

# Neuroticism score (20127)
y <- fread2(csv, select = "20127-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/neuroticism.rds")

# Forced expiratory volume in 1-second (FEV1), Best measure (20150)
y <- fread2(csv, select = "20150-0.0")[[1]]
y[log(y) < -1] <- NA
saveRDS(y, "data/ukbb-quant-pheno/FEV1.rds")

# Forced vital capacity (FVC), Best measure (20151)
y <- fread2(csv, select = "20151-0.0")[[1]]
y[log(y) < -0.5] <- NA
saveRDS(y, "data/ukbb-quant-pheno/FVC.rds")

# Fluid intelligence score (20191)
y <- fread2(csv, select = "20191-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/fluid_intelligence.rds")

# avMSE (20261)
y <- fread2(csv, select = "20261-0.0")[[1]]
y[abs(y) > 10] <- NA
saveRDS(y, "data/ukbb-quant-pheno/avMSE.rds")

# General happiness (20458)
y <- fread2(csv, select = "20458-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/less_happy.rds")

# General happiness with own health (20459)
y <- fread2(csv, select = "20459-0.0")[[1]]
y[y < 0] <- NA
saveRDS(y, "data/ukbb-quant-pheno/less_happy_with_health.rds")

# Body mass index (BMI) (21001)
y <- fread2(csv, select = "21001-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_BMI.rds")

# Mean carotid IMT (intima-medial thickness) at 120/150/210/240 degrees
y <- fread2(csv, select = c("22671-2.0", "22674-2.0",
                            "22677-2.0", "22680-2.0"))
y[y < 5.5] <- NA
y <- rowMeans(log(y))
saveRDS(y, "data/ukbb-quant-pheno/log_mean_carotid_IMT.rds")

# Body fat percentage (23099)
y <- fread2(csv, select = "23099-0.0")[[1]]
y <- y / 100
saveRDS(y, "data/ukbb-quant-pheno/fat_perc.rds")

# Whole body fat mass (23100)
y <- fread2(csv, select = "23100-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_fat_mass.rds")

# Whole body fat-free mass (23101)
y <- fread2(csv, select = "23101-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_fat_free_mass.rds")

# Whole body water mass (23102)
y <- fread2(csv, select = "23102-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_water_mass.rds")

# Impedance of whole body (23106)
y <- fread2(csv, select = "23106-0.0")[[1]]
y[log(y) < 5.8] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_impedance.rds")


#### ECG ####

# Ventricular rate (12336)
y <- fread2(csv, select = "12336-2.0")[[1]]
y[log(y) < 3 | log(y) > 5] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_ventricular_rate.rds")

# P duration (12338)
y <- fread2(csv, select = "12338-2.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/ECG_P_duration.rds")

# QRS duration (12340)
y <- fread2(csv, select = "12340-2.0")[[1]]
y[y < 40 | y > 150] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_ECG_QRS_duration.rds")

# PQ interval (22330)
y <- fread2(csv, select = "22330-2.0")[[1]]
y[log(y) > 6] <- NA
saveRDS(y, "data/ukbb-quant-pheno/ECG_PQ_interval.rds")

# QT interval (22331)
y <- fread2(csv, select = "22331-2.0")[[1]]
y[log(y) < 5.5] <- NA
saveRDS(y, "data/ukbb-quant-pheno/ECG_QT_interval.rds")

# QTC interval (22332)
y <- fread2(csv, select = "22332-2.0")[[1]]
y[y < 300 | y > 550] <- NA
saveRDS(y, "data/ukbb-quant-pheno/ECG_QTC_interval.rds")

# RR interval (22333)
y <- fread2(csv, select = "22333-2.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/ECG_RR_interval.rds")

# PP interval (22334)
y <- fread2(csv, select = "22334-2.0")[[1]]
y[log(y) < 6 | log(y) > 8] <- NA
saveRDS(y, "data/ukbb-quant-pheno/ECG_PP_interval.rds")


#### Blood count ####

# White blood cell (leukocyte) count (30000)
y <- fread2(csv, select = "30000-0.0")[[1]]
y[log(y) < 0 | log(y) > 4] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_leukocyte.rds")

# Red blood cell (erythrocyte) count (30010)
y <- fread2(csv, select = "30010-0.0")[[1]]
y[y < 2] <- NA
saveRDS(y, "data/ukbb-quant-pheno/erythrocyte.rds")

# Haemoglobin concentration (30020)
y <- fread2(csv, select = "30020-0.0")[[1]]
y[log(y) < 2] <- NA
saveRDS(y, "data/ukbb-quant-pheno/haemoglobin.rds")

# Haematocrit percentage (30030)
y <- fread2(csv, select = "30030-0.0")[[1]]
y <- y / 100
y[log(y) < -4] <- NA
saveRDS(y, "data/ukbb-quant-pheno/haematocrit_perc.rds")

# Mean corpuscular volume (30040)
y <- fread2(csv, select = "30040-0.0")[[1]]
y[y > 130] <- NA
saveRDS(y, "data/ukbb-quant-pheno/MCV.rds")

# Mean corpuscular haemoglobin (30050)
y <- fread2(csv, select = "30050-0.0")[[1]]
y[y < 10 | y > 50] <- NA
saveRDS(y, "data/ukbb-quant-pheno/MCH.rds")

# Red blood cell (erythrocyte) distribution width (30070)
y <- fread2(csv, select = "30070-0.0")[[1]]
y[y < 5 | y > 20] <- NA
saveRDS(y, "data/ukbb-quant-pheno/erythrocyte_width.rds")

# Platelet count (30080)
y <- fread2(csv, select = "30080-0.0")[[1]]
y[log(y) < 3] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_platelet.rds")

# Platelet crit (30090)
y <- fread2(csv, select = "30090-0.0")[[1]]
y[log(y) < -4] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_platelet_crit.rds")

# Mean platelet (thrombocyte) volume (30100)
y <- fread2(csv, select = "30100-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_platelet_volume.rds")

# Platelet distribution width (30110)
y <- fread2(csv, select = "30110-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_platelet_width.rds")

# Lymphocyte count (30120)
y <- fread2(csv, select = "30120-0.0")[[1]]
y[log(y) < -2 | log(y) > 3] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_lymphocyte.rds")

# Monocyte count (30130)
y <- fread2(csv, select = "30130-0.0")[[1]]
y[log(y) < -4 | log(y) > 2] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_monocyte.rds")

# Neutrophil count (30140)
y <- fread2(csv, select = "30140-0.0")[[1]]
y[log(y) < -1] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_neutrophil.rds")

# Lymphocyte percentage (30180)
y <- fread2(csv, select = "30180-0.0")[[1]]
y <- y / 100
y[log(y) < -4] <- NA
y[y > 0.7] <- NA
saveRDS(y, "data/ukbb-quant-pheno/lymphocyte_perc.rds")

# Monocyte percentage (30190)
y <- fread2(csv, select = "30190-0.0")[[1]]
y <- y / 100
y[log(y) < -4] <- NA
y[y > 0.2] <- NA
saveRDS(y, "data/ukbb-quant-pheno/monocyte_perc.rds")

# Neutrophil percentage (30200)
y <- fread2(csv, select = "30200-0.0")[[1]]
y <- y / 100
y[y < 0.25] <- NA
saveRDS(y, "data/ukbb-quant-pheno/neutrophil_perc.rds")

# Eosinophil percentage (30210)
y <- fread2(csv, select = "30210-0.0")[[1]]
y <- y / 100
y[log(y) < -4] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_eosinophil_perc.rds")

# Basophil percentage (30220)
y <- fread2(csv, select = "30220-0.0")[[1]]
y <- y / 100
y[log(y) < -4] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_basophil_perc.rds")

# Reticulocyte count (30250)
y <- fread2(csv, select = "30250-0.0")[[1]]
y[log(y) < -6 | log(y) > 0] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_reticulocyte.rds")

# Mean reticulocyte volume (30260)
y <- fread2(csv, select = "30260-0.0")[[1]]
y[y > 160] <- NA
saveRDS(y, "data/ukbb-quant-pheno/reticulocyte_volume.rds")

# Mean sphered cell volume (30270)
y <- fread2(csv, select = "30270-0.0")[[1]]
y[y < 60 | y > 120] <- NA
saveRDS(y, "data/ukbb-quant-pheno/sphered_cell_volume.rds")

# Immature reticulocyte fraction (30280)
y <- fread2(csv, select = "30280-0.0")[[1]]
y[y > 0.7] <- NA
saveRDS(y, "data/ukbb-quant-pheno/immature_reticulocyte_frac.rds")

# High light scatter reticulocyte count (30300)
y <- fread2(csv, select = "30300-0.0")[[1]]
y[log(y) < -6 | log(y) > -2] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_HLR_reticulocyte.rds")


#### Urine count ####

# Microalbumin in urine (30500)
y <- fread2(csv, select = "30500-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_microalbumin_urine.rds")

# Creatinine (enzymatic) in urine (30510)
y <- fread2(csv, select = "30510-0.0")[[1]]
y[log(y) < 6] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_creatinine_urine.rds")

# Potassium in urine (30520)
y <- fread2(csv, select = "30520-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_potassium_urine.rds")

# Sodium in urine (30530)
y <- fread2(csv, select = "30530-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/sodium_urine.rds")


#### Blood biochemistry ####

# Albumin (30600)
y <- fread2(csv, select = "30600-0.0")[[1]]
y[y < 30] <- NA
saveRDS(y, "data/ukbb-quant-pheno/albumin.rds")

# Alkaline phosphatase (30610)
y <- fread2(csv, select = "30610-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_ALP.rds")

# Alanine aminotransferase (30620)
y <- fread2(csv, select = "30620-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_ALT.rds")

# Apolipoprotein A (30630)
y <- fread2(csv, select = "30630-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/apoA.rds")

# Apolipoprotein B (30640)
y <- fread2(csv, select = "30640-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/apoB.rds")

# Aspartate aminotransferase (30650)
y <- fread2(csv, select = "30650-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_AST.rds")

# Urea (30670)
y <- fread2(csv, select = "30670-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_urea.rds")

# Calcium (30680)
y <- fread2(csv, select = "30680-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/calcium.rds")

# Cholesterol (30690)
y <- fread2(csv, select = "30690-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/cholesterol.rds")

# Creatinine (30700)
y <- fread2(csv, select = "30700-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_creatinine.rds")

# C-reactive protein (30710)
y <- fread2(csv, select = "30710-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_CRP.rds")

# Cystatin C (30720)
y <- fread2(csv, select = "30720-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_cystatinC.rds")

# Gamma glutamyltransferase (30730)
y <- log(fread2(csv, select = "30730-0.0")[[1]])
saveRDS(log(y), "data/ukbb-quant-pheno/log_log_gammaGT.rds")

# Glucose (30740)
y <- fread2(csv, select = "30740-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_glucose.rds")

# Glycated haemoglobin (HbA1c) (30750)
y <- fread2(csv, select = "30750-0.0")[[1]]
y[log(y) > 5] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_HbA1c.rds")

# HDL cholesterol (30760)
y <- fread2(csv, select = "30760-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_HDL.rds")

# IGF-1 (30770)
y <- fread2(csv, select = "30770-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_IGF1.rds")

# LDL direct (30780)
y <- fread2(csv, select = "30780-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/LDL.rds")

# Lipoprotein A (30790)
y <- fread2(csv, select = "30790-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_lipoA.rds")

# Oestradiol (30800)
y <- fread2(csv, select = "30800-0.0")[[1]]
y[sex == 1] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/F_log_oestradiol.rds")

y <- fread2(csv, select = "30800-0.0")[[1]]
y[sex == 0] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/M_log_oestradiol.rds")

# Phosphate (30810)
y <- fread2(csv, select = "30810-0.0")[[1]]
y[y > 2] <- NA
saveRDS(y, "data/ukbb-quant-pheno/phosphate.rds")

# Rheumatoid factor (30820)
y <- fread2(csv, select = "30820-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_rheumatoid_factor.rds")

# Sex hormone binding globulin (SHBG) (30830)
y <- fread2(csv, select = "30830-0.0")[[1]]
y[log(y) < 2] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/log_SHBG.rds")

# Total bilirubin (30840)
y <- fread2(csv, select = "30840-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_bilirubin.rds")

# Testosterone (30850)
y <- fread2(csv, select = "30850-0.0")[[1]]
y[sex == 1] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/F_log_testosterone.rds")

y <- fread2(csv, select = "30850-0.0")[[1]]
y[sex == 0] <- NA
saveRDS(log(y), "data/ukbb-quant-pheno/M_log_testosterone.rds")

# Total protein (30860)
y <- fread2(csv, select = "30860-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/protein.rds")

# Triglycerides (30870)
y <- fread2(csv, select = "30870-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_triglycerides.rds")

# Urate (30880)
y <- fread2(csv, select = "30880-0.0")[[1]]
saveRDS(y, "data/ukbb-quant-pheno/urate.rds")

# Vitamin D (30890)
y <- fread2(csv, select = "30890-0.0")[[1]]
saveRDS(log(y), "data/ukbb-quant-pheno/log_vitaminD.rds")


#### Also consider rank-based inverse normal transformation (RINT) ####

(pheno_files <- list.files("data/ukbb-quant-pheno", pattern = "^log_", full.names = TRUE))

rint <- function(x) qnorm(tail(dplyr::percent_rank(c(-Inf, Inf, x)), -2))

for (pheno_file in pheno_files) {
  saveRDS(rint(readRDS(pheno_file)), sub("log_", "rint_", pheno_file))
}

# and no transfo
for (pheno_file in pheno_files) {
  saveRDS(exp(readRDS(pheno_file)), sub("log_", "raw_", pheno_file))
}


#### Other binary UKBB phenotypes ####

bigassertr::assert_dir("data/ukbb-binary-pheno")

# Snoring (1210)
y <- fread2(csv, select = "1210-0.0")[[1]]
y <- ifelse(y < 0, NA, y == 1)
saveRDS(y, "data/ukbb-binary-pheno/snoring.rds")

# Daytime dozing / sleeping (narcolepsy) (1220)
y <- fread2(csv, select = "1220-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/narcolepsy.rds")

# Handedness (chirality/laterality) (1707)
y <- fread2(csv, select = "1707-0.0")[[1]]
y <- ifelse(y < 0 | y == 3, NA, y == 2)
saveRDS(y, "data/ukbb-binary-pheno/left_handed.rds")

# Hair colour (natural, before greying) (1747)
y <- fread2(csv, select = "1747-0.0")[[1]]
y <- ifelse(y < 0 | y == 6, NA, y == 2)
saveRDS(y, "data/ukbb-binary-pheno/red_hair.rds")

# Age first had sexual intercourse (2139)
y <- fread2(csv, select = "2139-0.0")[[1]]
y <- ifelse(y < 0 & y != -2, NA, y == -2)
saveRDS(y, "data/ukbb-binary-pheno/no_sex.rds")

# Wears glasses or contact lenses (2207)
y <- fread2(csv, select = "2207-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/glasses.rds")

# Plays computer games (2237)
y <- fread2(csv, select = "2237-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/play_computer.rds")

# Hearing difficulty/problems (2247)
y <- fread2(csv, select = "2247-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/bad_hearing.rds")

# Falls in the last year (2296)
y <- fread2(csv, select = "2296-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 1)
saveRDS(y, "data/ukbb-binary-pheno/fall_1y.rds")

# Fractured/broken bones in last 5 years (2463)
y <- fread2(csv, select = "2463-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/fracture_5y.rds")

# Length of menstrual cycle (3710)
y <- fread2(csv, select = "3710-0.0")[[1]]
y <- ifelse(y < 0 & y != -6, NA, y == -6)
y[sex == 1] <- NA
saveRDS(y, "data/ukbb-binary-pheno/F_irregular_menstrual_cycle.rds")

# Headaches for 3+ months (3799)
y <- fread2(csv, select = "3799-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/headaches_for_3m.rds")

# Ever depressed for a whole week (4598)
y <- fread2(csv, select = "4598-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/depressed_for_1w.rds")

# Alcohol drinker status (20117)
y <- fread2(csv, select = "20117-0.0")[[1]]
y <- ifelse(y < 0 | y == 1, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/drink_alcohol.rds")

# Ever smoked (20160)
y <- fread2(csv, select = "20160-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/ever_smoked.rds")

# Myopia diagnosis (20262)
y <- fread2(csv, select = "20262-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/myopia.rds")

# Ever addicted to any substance or behaviour (20401)
y <- fread2(csv, select = "20401-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/ever_addicted.rds")

# Depression possibly related to childbirth (20445)
y <- fread2(csv, select = "20445-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/postpartum_depression.rds")

# Depression possibly related to stressful or traumatic event (20447)
y <- fread2(csv, select = "20447-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/depression_stress_trauma.rds")

# Ever taken cannabis (20453)
y <- fread2(csv, select = "20453-0.0")[[1]]
y <- ifelse(y < 0 | y == 1, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/ever_cannabis.rds")


# Ever contemplated self-harm (20485)
# Recent thoughts of suicide or self-harm (20513)
y1 <- fread2(csv, select = "20485-0.0")[[1]]
y2 <- fread2(csv, select = "20513-0.0")[[1]]
y <- ifelse(y1 < 0 & y2 < 0, NA, y1 > 0 | y2 > 1)
saveRDS(y, "data/ukbb-binary-pheno/self_harm_thoughts.rds")

# Recent feelings of depression (20510)
y <- fread2(csv, select = "20510-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 1)
saveRDS(y, "data/ukbb-binary-pheno/depression_feelings.rds")

# Recent poor appetite or overeating (20511)
y <- fread2(csv, select = "20511-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 1)
saveRDS(y, "data/ukbb-binary-pheno/abnormal_appetite.rds")

# Recent feelings of foreboding (20512)
y <- fread2(csv, select = "20512-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 1)
saveRDS(y, "data/ukbb-binary-pheno/foreboding_feelings.rds")

# Trouble falling asleep (20533)
y <- fread2(csv, select = "20533-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/hard_falling_asleep.rds")

# Sleeping too much (20534)
y <- fread2(csv, select = "20534-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/sleeping_too_much.rds")

# Waking too early (20535)
y <- fread2(csv, select = "20535-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/waking_too_early.rds")

# Sensitive stomach (21064)
y <- fread2(csv, select = "21064-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/sensitive_stomach.rds")

# Diagnosed with coeliac disease or gluten sensitivity (21068)
y <- fread2(csv, select = "21068-0.0")[[1]]
y <- ifelse(y < 0, NA, y > 0)
saveRDS(y, "data/ukbb-binary-pheno/celiac_gluten.rds")
