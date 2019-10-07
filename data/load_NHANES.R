## Reading in the NHANES datasets

# Read in the Phtalate information
library(SASxport)
library(tidyverse)

setwd("~/SEM/data")

years = c("1516")

# functions to read the files
get_files <- function(str, years){
  files = system(paste("ls", str, sep=" "), intern=TRUE)
  tmp = files[map_dbl(years, function(y) ifelse(length(str_which(files, y)) == 1, str_which(files, y), NA))]
  tmp[!is.na(tmp)]
}

get_data <- function(files){
  map_dfr(files, function(f) as.data.frame(as.matrix(read.xport(f))))
}


# Get the data
# Possibility to add outcome variables related to dental status
strings = c("DEMO*.XPT","BMX*.XPT","BPX*.XPT","ALB_CR*.XPT",
            "APOB*.XPT","UADM_I*.XPT","UTAS_I*.XPT","CHLMDA_I*.XPT",
            "TCHOL_I*.XPT","CRCO_I*.XPT","CBC_I*.XPT","CUSEZN_I*.XPT",
            "COT_I*.XPT","UCOT_I*.XPT","DEET_I*.XPT","FERTIN_I",
            "FLDEP_I*.XPT","FLDEW_I*.XPT")

names = c("Demographics","Body_Measures","Blood_pressure",
          "Albumin_Creatinine_Urine","Apolipoprotein","Diamines",
          "Arsernic_Urine","Chlamydia","Cholesterol",
          "Chromium_Cobalt ","blood_count","copper_etc",
          "cotinine","cotinine_urine","deet","ferritin",
          "fluoride_plasma","fluoride_water")

data = tibble(names, strings) %>% 
  mutate(files = map(strings, function(s) get_files(s, years)),
         data = map(files, get_data))

# Demographics
# https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Demographics&CycleBeginYear=2015
people_nhanes = data %>%
  filter(names == "Demographics") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(c(RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC3, DMDEDUC2,
           DMDMARTL,INDFMPIR,RIDEXPRG)) 
# RIDAGEYR: Age in years at screening
# RIAGENDR: Gender
# RIDRETH1: Race/Hispanic origin 
# INDFMPIR: Ratio of family income to poverty
# DMDEDUC3: Education level - Children/Youth 6-19
# DMDEDUC2: Education level - Adults 20+
# DMDMARTL: Marital status
# RIDEXPRG: Pregnancy status at exam

# For chemicals
#https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Laboratory&CycleBeginYear=2015

# Fluoride - Water (FLDEW_I)
fluoride_water  = data %>%
  filter(names == "fluoride_plasma") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBDWFL)
# LBDWFL - Fluoride, water (mg/L) average 2 values

# Fluoride - Plasma (FLDEP_I)
fluoride_plasma  = data %>%
  filter(names == "fluoride_plasma") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBDPFL)
# LBDPFL - Fluoride, plasma (umol/L) average 2

# Ferritin (FERTIN_I)
ferritin = data %>%
  filter(names == "ferritin") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXFER)
# LBXFER - Ferritin (ng/mL)


# DEET Metabolite - Urine (DEET_I)
deet = data %>%
  filter(names == "deet") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXDEA)
# URXDEA: DEET acid (ng/mL)

# Cotinine, Hydroxycotinine, & Other Nicotine Metabolites and Analogs - Urine (UCOT_I)
cotinine_urine = data %>%
  filter(names == "cotinine_urine") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXCOTT,URXHCTT,URXANBT,URXANTT,
         URXCOXT,URXHPBT,URXNICT,URXNNCT,
         URXNOXT,URDTNE2,URDTNE3,URDTNE6,URDTNE7)
# URXCOTT: Total Cotinine, urine (ng/mL)
# URXHCTT: Total Hydroxycotinine, urine (ng/mL)
# URXANBT: Anabasine, urine (ng/mL)
# URXANTT: Anatabine, urine (ng/mL)
# URXCOXT: Cotinine-n-oxide, urine (ng/mL)
# URXHPBT: 1-(3P)-1-but-4-carbox acid (ng/mL)
# URXNICT: Nicotine, urine (ng/mL)
# URXNNCT: Nornicotine, urine (ng/mL)
# URXNOXT: Nicotine-1 N-oxide, urine (ng/mL)
# URDTNE2: TNE - 2 (nmol/mL)
# URDTNE3: TNE - 3 (nmol/mL)
# URDTNE6: TNE - 3 (nmol/mL)
# URDTNE7: TNE - 3 (nmol/mL)


# Cotinine and Hydroxycotinine - Serum (COT_I)
cotinine = data %>%
  filter(names == "cotinine") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXCOT,LBXHCT)
# LBXCOT: Cotinine, Serum (ng/mL)
# LBXHCT: Hydroxycotinine, Serum (ng/mL)


# Copper, Selenium & Zinc - Serum (CUSEZN_I)
copper_etc = data %>%
  filter(names == "copper_etc") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXBCO)
# LBXBCR: Chromium (ug/L)
# LBXBCO: Cobalt (ug/L)

# Chromium & Cobalt (CRCO_I)
chromium_cobalt = data %>%
  filter(names == "Chromium_Cobalt") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXSCU,LBXSSE,LBXSZN)
# LBXSCU: Serum Copper (ug/dL)
# LBXSSE: Serum Selenium (ug/L)
# LBXSZN: Serum Zinc (ug/dL)


# Cholesterol - Total (TCHOL_I)
chol = data %>%
  filter(names == "Cholesterol") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXTC)
# LBXTC: Total Cholesterol (mg/dL)

# Arsenic - Total - Urine (UTAS_I)
arsenic = data %>%
  filter(names == "Arsernic_Urine") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXUAS)
# URXUAS: Urinary arsenic, total (ug/L)

#Chlamydia - Urine (CHLMDA_I)
chlamydia = data %>%
  filter(names == "Chlamydia") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXUCL)
# URXUCL - Chlamydia, Urine

# Albumin & Creatinine - Urine
albumin_creatinine_nhanes = data %>%
  filter(names == "Albumin_Creatinine_Urine") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXUCR,URXUMA)
# URXUMS: Albumin, urine (mg/L)
# URXUCR: Creatinine, urine (mg/dL)

#Apolipoprotein
apolipoprotein_nhanes = data %>%
  filter(names == "Apolipoprotein") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXAPB)
# LBXAPB: Apolipoprotein (B) (mg/dL)


# Aromatic Diamines - Urine (UADM_I)
diamines_nhanes = data %>%
  filter(names == "Diamines") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URX4TDA,URX6TDA,URX4MDA,URX5NDA,
         URXPPDA)
# URX4TDA: 2,4-Diaminotoluene (4TDA) (ng/mL)
# URX6TDA: 2,6-Diaminotoluene (6TDA) (ng/mL)
# URX4MDA: 4MDA (ng/mL)
# URX5NDA: 1,5-Diaminonaphthalene (5NDA) (ng/mL)
# URXPPDA: p-Phenylenediamine (PPDA) (ng/mL)







# Potential Outcomes
# Complete Blood Count with 5-Part Differential - Whole Blood (CBC_I)
blood_count = data %>%
  filter(names == "blood_count") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXWBCSI,LBXLYPCT,LBXMOPCT,LBXNEPCT,
         LBXEOPCT,LBXBAPCT,LBXRBCSI,LBXHGB,
         LBXHCT,LBXMCVSI,LBXRDW,LBXMPSI)
# LBXWBCSI: White blood cell count (1000 cells/uL)
# LBXLYPCT: Lymphocyte percent (%)
# LBXMOPCT: Monocyte percent (%)
# LBXNEPCT: Segmented neutrophils percent (%)
# LBXEOPCT: Eosinophils percent (%)
# LBXBAPCT: Basophils percent (%)
# LBXRBCSI: Red blood cell count (million cells/uL)
# LBXHGB: Hemoglobin (g/dL)
# LBXHCT: Hematocrit (%)
# LBXMCVSI: Mean cell volume (fL)
# LBXRDW: Red cell distribution width (%)
# LBXMPSI: Mean platelet volume (fL)

























phthalates_nhanes = c("URXMBP", "URXMIB", "URXMEP", 
                      "URXMZP", "URXMCP", "URXECP", 
                      "URXMHH", "URXMOH", "URXMHP")
# URXMBP: Mono-n-butyl phthalate (ng/mL)
# URXMIB: Mono-isobutyl phthalate
# URXMEP: Mono-ethyl phthalate (ng/mL)
# URXMZP: Mono-benzyl phthalate (ng/mL)
# URXMCP: Mono-cyclohexyl phthalate (ng/mL)
# URXECP: Mono-2-ethyl-5-carboxypentyl phthalate
# URXMHH: Mono-(2-ethyl-5-hydroxyhexyl) phthalate
# URXMOH: Mono-(2-ethyl-5-oxohexyl) phthalate
# URXMHP: Mono-(2-ethyl)-hexyl phthalate (ng/mL)


chem_nhanes = data %>%
  filter(names == "Phthlates") %>%
  select(data) %>%
  unnest() %>%
  select(SEQN, phthalates_nhanes) %>%
  filter_all(all_vars(!is.na(.)))
# SEQN: Respondent sequence number


bmi_nhanes = data %>%
  filter(names == "Body_Measures") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(BMXBMI, BMXWAIST)
# BMXBMI: Body Mass Index (kg/m**2)
# BMXWAIST: Waist Circumference (cm)

# Only 2009-2010 and 2011-2012
pesticides_nhanes = data %>%
  filter(names == "Organochlorine_Pesticides") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URX14D,URXDCB)
# URX14D: 2,5-dichlorophenol (ug/L) result
# URXDCB: 2,4-dichlorophenol (ug/L) result


# people that have Phalates do not have metals
metals_nhanes = data %>%
  filter(names == "Metals") %>%
  select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXUBA,URXUBE,URXUCD,
         URXUCO,URXUCS,URXUMO,
         URXUPB,URXUPT,URXUSB,
         URXUTL,URXUTU,URXUUR)
# URXUBA: Barium, urine (ug/L)
# URXUBE: Beryllium, urine (ug/L)
# URXUCD: Cadmium, urine (ug/L)
# URXUCO: Cobalt, urine (ug/L)
# URXUCS: Cesium, urine (ug/L)
# URXUMO: Molybdenum, urine (ug/L)
# URXUPB: Lead, urine (ug/L)
# URXUPT: Platinum, urine (ug/L)
# URXUSB: Antimony, urine (ug/L)
# URXUTL: Thallium, urine (ug/L)
# URXUTU: Tungsten, urine (ug/L)
# URXUUR: Uranium, urinary (ug/L)

# Missing year 2013-2014
phenols_nhanes = data %>%
  filter(names == "Phenols") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(URXBPH,URX4TO,URXTRS,URXBP3,
         URXPPB,URXBUP,URXEPB,URXMPB)
# URXBPH: Urinary Bisphenol A
# URX4TO: Urinary 4-tert-octyl 
# URXTRS: Urinary 2,4,4’-Trichloro-2’-hydroxyphenyl ether (Triclosan)		
# URXBP3: Urinary 2-Hydroxy-4-metoxybenzophenone (Benzophenone-3)	
# URXPPB: Urinary Propyl paraben	
# URXBUP: Urinary Butyl paraben		
# URXEPB: Urinary Ethyl paraben	
# URXMPB: Urinary Methyl paraben	

chol_nhanes = data %>%
  filter(names == "Cholesterol") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(LBXTC)
# LBXTC: Total Cholesterol (mg/dL)


data_nhanes = data.frame(chem_nhanes,people_nhanes,bmi_nhanes,
                         creatinine_nhanes,chol_nhanes,pesticides_nhanes)

# Rename the variables=
colnames(data_nhanes)[colnames(data_nhanes)=="URXUCR"] = "Creatinine"
colnames(data_nhanes)[colnames(data_nhanes)=="LBXTC"] = "TotChol"
colnames(data_nhanes)[colnames(data_nhanes)=="RIDAGEYR"] = "Age"
colnames(data_nhanes)[colnames(data_nhanes)=="RIAGENDR"] = "Gender"
colnames(data_nhanes)[colnames(data_nhanes)=="SEQN"] = "ID"
colnames(data_nhanes)[colnames(data_nhanes)=="BMXBMI"] = "BMI"
colnames(data_nhanes)[colnames(data_nhanes)=="BMXWAIST"] = "WAIST"
colnames(data_nhanes)[colnames(data_nhanes)=="INDFMPIR"] = "Ratio_income_poverty"
colnames(data_nhanes)[colnames(data_nhanes)=="RIDRETH1"] = "Race"

quantile(data_nhanes$URXMBP)
quantile(data_nhanes$URXMIB)

#log transform
data_nhanes$URXMBP = log10(data_nhanes$URXMBP)
data_nhanes$URXMIB = log10(data_nhanes$URXMIB)
data_nhanes$URXMEP = log10(data_nhanes$URXMEP)
data_nhanes$URXMZP = log10(data_nhanes$URXMZP)
data_nhanes$URXMCP = log10(data_nhanes$URXMCP)
data_nhanes$URXECP = log10(data_nhanes$URXECP)
data_nhanes$URXMHH = log10(data_nhanes$URXMHH)
data_nhanes$URXMOH = log10(data_nhanes$URXMOH)
data_nhanes$URXMHP = log10(data_nhanes$URXMHP)
data_nhanes$URX14D = log10(data_nhanes$URX14D)
data_nhanes$URXDCB = log10(data_nhanes$URXDCB)
data_nhanes$Creatinine = log10(data_nhanes$Creatinine)
data_nhanes$TotChol = log10(data_nhanes$TotChol)

# save data
save(data_nhanes, file = "nhanes_complete.RData")

#log transform + creatanine adjusted
log_creat = log(data_nhanes$Creatinine)
data_nhanes$URXMBP = data_nhanes$URXMBP - log_creat
data_nhanes$URXMIB = data_nhanes$URXMIB - log_creat
data_nhanes$URXMEP = data_nhanes$URXMEP - log_creat
data_nhanes$URXMZP = data_nhanes$URXMZP - log_creat
data_nhanes$URXMCP = data_nhanes$URXMCP - log_creat
data_nhanes$URXECP = data_nhanes$URXECP - log_creat
data_nhanes$URXMHH = data_nhanes$URXMHH - log_creat
data_nhanes$URXMOH = data_nhanes$URXMOH - log_creat
data_nhanes$URXMHP = data_nhanes$URXMHP - log_creat
data_nhanes$URX14D = data_nhanes$URX14D - log_creat
data_nhanes$URXDCB = data_nhanes$URXDCB - log_creat

# save log transformed
save(data_nhanes, file = "nhanes_complete_adj.RData")




