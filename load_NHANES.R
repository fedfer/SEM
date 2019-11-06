## Reading in the NHANES datasets

# Missing body info

library(SASxport)
library(tidyverse)
library(plyr)

setwd("~/SEM/data_sas")

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
            "COT_I*.XPT","UCOT_I*.XPT","DEET_I*.XPT","FERTIN_I*.XPT",
            "FLDEP_I*.XPT","FLDEW_I*.XPT","FOLATE_I*.XPT","FOLFMS_I*.XPT",
            "GHB_I*.XPT","SSNEON_I*.XPT","UIO_I*.XPT",
            "PBCD_I*.XPT","UHG_I*.XPT","IHGEM_I*.XPT","UM_I*.XPT",
            "PFAS_I*.XPT","PHTHTE_I*.XPT","TST_I*.XPT","UAS_I*.XPT",
            "BIOPRO_I*.XPT","TRICH_I*.XPT","UVOC_I*.XPT",
            "VOCWB_I*.XPT")


names = c("demographics","body_measures","blood_pressure",
          "albumin_creatinine_urine","apolipoprotein","diamines",
          "arsernic_urine","chlamydia","cholesterol",
          "chromium_cobalt","blood_count","copper_etc",
          "cotinine","cotinine_urine","deet","ferritin",
          "fluoride_plasma","fluoride_water","folate","folate_serum",
          "glycohemoglobin","neonicotinoids","iodine","metals_blood",
          "mercury_urine","mercury_blood","metals_urine","pfas",
          "phthalate","steroid","spec_arsernic_urine","bio",
          "trichomonas","voc_urine","voc_blood")

data = tibble(names, strings) %>% 
  mutate(files = map(strings, function(s) get_files(s, years)),
         data = map(files, get_data))

#### DEMOGRAPHICS ####
# https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Demographics&CycleBeginYear=2015
people = data %>%
  filter(names == "demographics") %>%
  dplyr::select(data) %>%
  unnest() %>%
  dplyr::select(SEQN,RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC3, DMDEDUC2,
           DMDMARTL,INDFMPIR,RIDEXPRG)
# RIDAGEYR: Age in years at screening
# RIAGENDR: Gender
# RIDRETH1: Race/Hispanic origin 
# INDFMPIR: Ratio of family income to poverty
# DMDEDUC3: Education level - Children/Youth 6-19
# DMDEDUC2: Education level - Adults 20+
# DMDMARTL: Marital status
# RIDEXPRG: Pregnancy status at exam

# For chemicals and effect modifiers
#https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Laboratory&CycleBeginYear=2015


#### CHEMICALS ####
# Volatile Organic Compounds and Trihalomethanes/MTBE - Blood (VOCWB_I)
voc_blood = data %>%
  filter(names == "voc_blood") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBX2DF,LBX4CE,LBXV06,LBXV07N,LBXV08N,
         LBXV3B,LBXV4C,LBXVBF,LBXVBM,LBXVBZ,
         LBXVBZN,LBXVC6,LBXVCB,LBXVCF,LBXVCM,LBXVCT,
         LBXVDB,LBXVDE,LBXVDEE,LBXVDX,LBXVEA,LBXVEB,LBXVEC,
         LBXVFN,LBXVIBN,LBXVIPB,LBXVMC,LBXVMCP,LBXVME,LBXVNB,
         LBXVOX,LBXVTC,LBXVTE,LBXVTFT,LBXVTHF,LBXVTO,LBXVTP,
         LBXVVB,LBXVXY
         )
# LBX2DF Blood 2,5-Dimethylfuran (ng/mL)
# LBX4CE: Blood 1,1,1,2-Tetrachloroethane (ng/mL)
# LBXV06: Blood Hexane (ng/mL)
# LBXV07N: Blood Heptane (ng/mL)
# LBXV08N: Blood Octane (ng/mL)
# LBXV2A: Blood 1,2-Dichloroethane (ng/mL)
# LBXV3B: Blood 1,3-Dichlorobenzene (ng/mL)
# LBXV4C: Blood Tetrachloroethene (ng/mL)
# LBXVBF: Blood Bromoform (ng/mL)
# LBXVBM: Blood Bromodichloromethane (ng/mL)
# LBXVBZ: Blood Benzene (ng/mL)
# LBXVBZN: Blood Benzonitrile (ng/mL) 
# LBXVC6: Blood Cyclohexane (ng/mL)
# LBXVCB: Blood Chlorobenzene (ng/mL)
# LBXVCF: Blood Chloroform (ng/mL)
# LBXVCM: Blood Dibromochloromethane (ng/mL)
# LBXVCT: Blood Carbon Tetrachloride (ng/mL)
# LBXVDB: Blood 1,4-Dichlorobenzene (ng/mL)
# LBXVDE: Blood 1,2-Dibromoethane (ng/ml)
# LBXVDEE: Blood Diethyl Ether (ng/mL)
# LBXVDX: Blood 1,4-Dioxane (ng/mL)
# LBXVEA: Blood Ethyl Acetate (ng/mL)
# LBXVEB: Blood Ethylbenzene (ng/mL)
# LBXVEC: Blood Chloroethane (ng/mL)
# LBXVFN: Blood Furan (ng/ml)
# LBXVIBN: Blood Isobutyronitrile (ng/mL)
# LBXVIPB: Blood Isopropylbenzene (ng/ml)
# LBXVMC: Blood Methylene Chloride (ng/mL)
# LBXVMCP: Blood Methylcyclopentane (ng/mL)
# LBXVME: Blood MTBE (ng/mL)
# LBXVNB: Blood Nitrobenzene (ng/mL)
# LBXVOX: Blood o-Xylene (ng/mL)
# LBXVTC: Blood Trichloroethene (ng/mL)
# LBXVTE: Blood 1,1,1-Trichloroethane (ng/mL)
# LBXVTFT: Blood ααα-Trifluorotoluene (ng/mL)
# LBXVTHF: Blood Tetrahydrofuran (ng/mL)
# LBXVTO: Blood Toluene (ng/mL)
# LBXVTP: Blood 1,2,3-Trichloropropane (ng/ml)
# LBXVVB: Blood Vinyl Bromide (ng/mL)
# LBXVXY: Blood m-/p-Xylene (ng/mL)

# Volatile Organic Compound (VOC) Metabolites - Urine (UVOC_I)
voc_urine = data %>%
  filter(names == "voc_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXAAM,URXAMC,URXATC,URXBMA,URXBPM,
         URXCEM,URXCYHA,URXCYM,URX1DC,
         URX2DC,URXDHB,URXDPM,URXGAM,
         URXHPM,URXHP2,URXIPM1,URXIPM3,URXMAD,
         URX2MH,URX34M,URXMB1,URXMB2,URXMB3,
         URXPHE,URXPHG,URXPMA,URXPMM,URXTCV)
# URXAAM: N-Acetyl-S-(2-carbamoylethyl)-L-cysteine (ng/mL)
# URXAMC: N-Acetyl-S-(N-methylcarbamoyl)-L-cysteine (ng/mL)
# URXATC: 2-Aminothiazoline-4-carboxylic acid (ng/mL)
# URXBMA: N-Acetyl-S-(benzyl)-L-cysteine (ng/mL)
# URXBPM: N-Acetyl-S-(n-propyl)-L-cysteine (ng/mL)
# URXCEM: N-Acetyl-S-(2-carboxyethyl)-L-cysteine (ng/mL)
# URXCYHA: N-Acetyl-S-(1-cyano-2-hydroxyethyl)-L-cysteine (ng/mL)
# URXCYM: N-Acetyl-S-(2-cyanoethyl)-L-cysteine (ng/mL)
# URX1DC: N-Acetyl-S-(1,2-dichlorovinyl)-L-cysteine (ng/mL)
# URX2DC: N-Acetyl-S-(2,2-dichlorovinyl)-L-cysteine (ng/mL)
# URXDHB: N-Acetyl-S-(3,4-dihydroxybutyl)-L-cysteine (ng/mL)
# URXDPM: N-Acetyl-S-(dimethylphenyl)-L-cysteine (ng/mL)
# URXGAM: N-Acetyl-S-(2-carbamoyl-2-hydroxyethyl)-L-cysteine (ng/mL)
# URXHPM: N-Acetyl-S-(3-hydroxypropyl)-L-cysteine (ng/mL)
# URXHP2: N-Acetyl-S-(2-hydroxypropyl)-L-cysteine (ng/mL)
# URXIPM1: N-Acetyl-S-(2-hydroxy-3-methyl-3-butenyl)-L-cysteine +N-Acetyl-S-(2-hydroxy-2-methyl-3-butenyl)-L-cysteine
# URXIPM3: N-Acetyl- S- (4- hydroxy- 2- methyl- 2- butenyl)- L- cysteine (ng/mL)
# URXMAD Mandelic acid (ng/mL)
# URX2MH: 2-Methylhippuric acid (ng/mL)
# URX34M: 3- and 4-Methylhippuric acid (ng/mL)
# URXMB1: N-Acetyl-S-(1-hydroxymethyl-2-propenyl)-L-cysteine (ng/mL)
# URXMB2: N-Acetyl-S-(2-hydroxy-3-butenyl)-L-cysteine (ng/mL)
# URXMB3: N-Acetyl-S-(4-hydroxy-2-butenyl)-L-cysteine (ng/mL)
# URXPHE: N-Acetyl-S-(phenyl-2-hydroxyethyl)-L-cysteine (ng/mL)
# URXPHG: Phenylglyoxylic acid (ng/mL)
# URXPMA: N-Acetyl-S-(phenyl)-L-cysteine (ng/mL)
# URXPMM: N-Acetyl-S-(3-hydroxypropyl-1-methyl)-L-cysteine (ng/mL)
# URXTCV: N-Acetyl-S-(trichlorovinyl)-L-cysteine (ng/mL)

# Trichomonas - Urine (TRICH_I)
trichomonas = data %>%
  filter(names == "trichomonas") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUTRI) 
# URXUTRI: Trichomonas, Urine

# Standard Biochemistry Profile (BIOPRO_I)
bio = data %>%
  filter(names == "bio") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXSATSI,LBXSAL,LBXSAPSI,LBXSASSI,
         LBXSC3SI,LBXSBU,LBXSCLSI,LBXSCH,
         LBXSCK,LBXSCR,LBXSGTSI,LBXSGL,LBXSIR,
         LBXSPH,LBXSKSI,LBXSNASI,LBXSTB,LBXSCA,LBXSTP,
         LBXSTR,LBXSUA) 
# LBXSATSI: Alanine Aminotransferase (ALT) (U/L)
# LBXSAL: Albumin, refrigerated serum (g/dL)
# LBXSAPSI: Alkaline Phosphatase (ALP) (U/L)
# LBXSASSI: Aspartate Aminotransferase (AST) (U/L)
# LBXSC3SI: Bicarbonate (mmol/L)
# LBXSBU: Blood Urea Nitrogen (g/dL)
# LBXSCLSI: Chloride (mmol/L)
# LBXSCH: Cholesterol, refrigerated serum (mg/dL)
# LBXSCK: Creatine Phosphokinase (CPK) (IU/L)
# LBXSCR: Creatinine, refrigerated serum (mg/dL)
# LBXSGTSI: Gamma Glutamyl Transferase (GGT) (U/L)
# LBXSGL: Glucose, refrigerated serum (mg/dL)
# LBXSIR: Iron, refrigerated serum (ug/dL)
# LBXSPH: Phosphorus (mg/dL)
# LBXSKSI: Potassium (mmol/L)
# LBXSNASI: Sodium (mmol/L)
# LBXSTB: Total Bilirubin (mg/dL)
# LBXSCA: Total Calcium (mg/dL)
# LBXSTP: Total Protein (g/dL)
# LBXSTR: Triglycerides, refrig serum (mg/dL)
# LBXSUA: Uric Acid (mg/dL)

# Speciated Arsenics - Urine (UAS_I)
spec_arsernic_urine = data %>%
  filter(names == "spec_arsernic_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  dplyr::select(SEQN,URXUAS3,URXUAS5,URXUAB,URXUAC,
         URXUDMA,URXUMMA)
# URXUAS3: Urinary Arsenous Acid
# URXUAS5: Urinary Arsenic acid
# URXUAB: Urinary Arsenobetaine
# URXUAC: Urinary Arsenocholine
# URXUDMA: Urinary Dimethylarsinic Acid
# URXUMMA: Urinary Monomethylarsonic Acid


# Sex Steroid Hormone - Serum (TST_I)
steroid  = data %>%
  filter(names == "steroid") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXTST,LBXEST,LBXSHBG)
# LBXTST: Testosterone, total (ng/dL)
# LBXEST: Estradiol (pg/mL)
# LBXSHBG: SHBG (nmol/L)

# Phthalates and Plasticizers Metabolites - Urine (PHTHTE_I)
phthalate  = data %>%
  filter(names == "phthalate") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXCNP,URXCOP,URXECP,URXHIBP,URXMBP,URXMC1,
         URXMCOH,URXMEP,URXMHBP,URXMHH,URXMHNC,URXMHP,
         URXMIB,URXMNP,URXMOH,URXMZP)
# URXCNP:	Mono(carboxyisononyl) phthalate (ng/mL)	
# URXCOP:	Mono(carboxyisoctyl) phthalate (ng/mL)
# URXECP:	Mono-2-ethyl-5-carboxypentyl phthalate (ng/mL)	
# URXHIBP:	MHIBP phthalate (ng/mL)	
# URXMBP:	Mono-n-butyl phthalate (ng/mL)	
# URXMC1:	Mono-(3-carboxypropyl) phthalate (ng/mL)	
# URXMCOH:	MCOCH phthalate (ng/mL)	
# URXMEP:	Mono-ethyl phthalate (ng/mL)	
# URXMHBP:	Mono-3-hydroxy-n-butyl phthalate (ng/mL)	
# URXMHH:	Mono-(2-ethyl-5-hydroxyhexyl) phthalate (ng/mL)	
# URXMHNC:	Cyclohexane 1,2-dicarboxylic acid monohydroxy isononyl ester (ng/mL)	
# URXMHP:	Mono-(2-ethyl)-hexyl phthalate (ng/mL)	
# URXMIB:	Mono-isobutyl phthalate (ng/mL)	
# URXMNP:	Mono-isononyl phthalate (ng/mL)	
# URXMOH:	Mono-(2-ethyl-5-oxohexyl) phthalate (ng/mL)	
# URXMZP:	Mono-benzyl phthalate (ng/mL)

# Perfluoroalkyl and Polyfluoroalkyl (PFAS_I)
pfas  = data %>%
  filter(names == "pfas") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXPFDE,LBXPFHS,LBXMPAH,LBXPFNA,LBXPFUA,
         LBXPFDO,LBXNFOA,LBXBFOA,LBXNFOS,LBXMFOS)
# LBXPFDE: Perfluorodecanoic acid (ng/mL)
# LBXPFHS: Perfluorohexane sulfonic acid (ng/mL)
# LBXMPAH: 2-(N-methyl-PFOSA)acetic acid (ng/mL)
# LBXPFNA: Perfluorononanoic acid (ng/mL)
# LBXPFUA: Perfluoroundecanoic acid (ng/mL)
# LBXPFDO: Perfluorododecanoic acid (ng/mL)
# LBXNFOA: n-perfluorooctanoic acid (ng/mL)
# LBXBFOA: Br. perfluorooctanoic acid iso (ng/mL)
# LBXNFOS: n-perfluorooctane sulfonic acid (ng/mL)
# LBXMFOS: Sm-PFOS (ng/mL)

# Metals - Urine (UM_I)
metals_urine  = data %>%
  filter(names == "metals_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUBA,URXUCD,URXUCO,URXUCS,URXUMO,
         URXUMN,URXUPB,URXUSB,URXUSN,URXUSR,
         URXUTL,URXUTU,URXUUR)
# URXUBA: Barium, urine (ug/L)
# URXUCD: Cadmium, urine (ug/L)
# URXUCO: Cobalt, urine (ug/L)
# URXUCS: Cesium, urine (ug/L)
# URXUMO: Molybdenum, urine (ug/L)
# URXUMN: Manganese, urine (ug/L)
# URXUPB: Lead, urine (ug/L)
# URXUSB: Antimony, urine (ug/L)
# URXUSN: Tin, urine (ug/L)
# URXUSR: Strontium, urine (ug/L)
# URXUTL: Thallium, urine (ug/L)
# URXUTU: Tungsten, urine (ug/L)
# URXUUR: Uranium, urine (ug/L)

# Mercury: Inorganic, Ethyl and Methyl – Blood (IHGEM_I)
mercury_blood  = data %>%
  filter(names == "mercury_blood") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXIHG,LBXBGE,LBXBGM)
# LBXIHG: Mercury, inorganic (ug/L)
# LBXBGE: Mercury, ethyl (ug/L)
# LBXBGM: Mercury, methyl (ug/L)


# Mercury - Urine (UHG_I)
mercury_urine  = data %>%
  filter(names == "mercury_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUHG)
# URXUHG: Urine Mercury (ng/mL)

# Lead, Cadmium, Total Mercury, Selenium & Manganese - Blood (PBCD_I)
metals_blood  = data %>%
  filter(names == "metals_blood") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXBPB,LBXBCD,LBXTHG,LBXBSE,LBXBMN)
# LBXBPB: Blood lead (ug/dL)
# LBXBCD: Blood cadmium (ug/L)
# LBXTHG: Blood mercury, total (ug/L)
# LBXBSE: Blood selenium (ug/L)
# LBXBMN: Blood manganese (ug/L)

# Iodine - Urine (UIO_I)
iodine  = data %>%
  filter(names == "iodine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUIO)
# URXUIO: Iodine, urine (ng/mL)

# Neonicotinoids - Urine - Surplus (SSNEON_I)
neonicotinoids  = data %>%
  filter(names == "neonicotinoids") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,SSIMID,SSACET,SSCLOT,
         SSTHIA,SSOHIM,SSAND)
# SSIMID: Imidacloprid (ug/L)
# SSACET: Acetamiprid (ug/L)
# SSCLOT: Clothianidin (ug/L)
# SSTHIA: Thiacloprid (ug/L)
# SSOHIM: 5-Hydroxyimidacloprid (ug/L)
# SSAND - N-Desmethylacetamiprid (ug/L)

# Glycohemoglobin (GHB_I)
glycohemoglobin  = data %>%
  filter(names == "glycohemoglobin") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXGH)
# LBXGH: Glycohemoglobin (%)

#Folate Forms - Total & Individual - Serum (FOLFMS_I)
folate_serum  = data %>%
  filter(names == "folate_serum") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBDFOT,LBDSF1LC,LBDSF2LC,LBDSF3LC,
         LBDSF4LC,LBDSF5LC,LBDSF6LC)
# LBDFOT: Serum total folate (ng/mL)
# LBDSF1LC: 5-Methyl-tetrahydrofolate cmt
# LBDSF2LC: Folic acid cmt
# LBDSF3LC: 5-Formyl-tetrahydrofolate cmt
# LBDSF4LC: Tetrahydrofolate cmt
# LBDSF5LC: 5,10-Methenyl-tetrahydrofolate cmt
# LBDSF6LC: Mefox oxidation product cmt

# Folate - RBC (FOLATE_I)
folate  = data %>%
  filter(names == "folate") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBDRFO)
# LBDRFO: RBC folate (ng/mL)

# Fluoride - Water (FLDEW_I)
fluoride_water  = data %>%
  filter(names == "fluoride_water") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBDWFL)
# LBDWFL - Fluoride, water (mg/L) average 2 values

# Fluoride - Plasma (FLDEP_I)
fluoride_plasma  = data %>%
  filter(names == "fluoride_plasma") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBDPFL)
# LBDPFL - Fluoride, plasma (umol/L) average 2

# Ferritin (FERTIN_I)
ferritin = data %>%
  filter(names == "ferritin") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXFER)
# LBXFER - Ferritin (ng/mL)


# DEET Metabolite - Urine (DEET_I)
deet = data %>%
  filter(names == "deet") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXDEA)
# URXDEA: DEET acid (ng/mL)

# Cotinine, Hydroxycotinine, & Other Nicotine Metabolites and Analogs - Urine (UCOT_I)
cotinine_urine = data %>%
  filter(names == "cotinine_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXCOTT,URXHCTT,URXANBT,URXANTT,
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
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXCOT,LBXHCT)
# LBXCOT: Cotinine, Serum (ng/mL)
# LBXHCT: Hydroxycotinine, Serum (ng/mL)

# Copper, Selenium & Zinc - Serum (CUSEZN_I)
copper_etc = data %>%
  filter(names == "copper_etc") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXSCU,LBXSSE,LBXSZN)
# LBXSCU: Serum Copper (ug/dL)
# LBXSSE: Serum Selenium (ug/L)
# LBXSZN: Serum Zinc (ug/dL)
  
# Chromium & Cobalt (CRCO_I)
chromium_cobalt = data %>%
  filter(names == "chromium_cobalt") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXBCR,LBXBCO)
# LBXBCR: Chromium (ug/L)
# LBXBCO: Cobalt (ug/L)


# Cholesterol - Total (TCHOL_I)
chol = data %>%
  filter(names == "cholesterol") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXTC)
# LBXTC: Total Cholesterol (mg/dL)

# Arsenic - Total - Urine (UTAS_I)
arsenic = data %>%
  filter(names == "arsernic_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUAS)
# URXUAS: Urinary arsenic, total (ug/L)

#Chlamydia - Urine (CHLMDA_I)
chlamydia = data %>%
  filter(names == "chlamydia") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUCL)
# URXUCL - Chlamydia, Urine

# Albumin & Creatinine - Urine
albumin_creatinine_nhanes = data %>%
  filter(names == "albumin_creatinine_urine") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URXUMA,URXUCR)
# URXUMA: Albumin, urine (mg/L)
# URXUCR: Creatinine, urine (mg/dL)

#Apolipoprotein
apolipoprotein = data %>%
  filter(names == "apolipoprotein") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXAPB)
# LBXAPB: Apolipoprotein (B) (mg/dL)


# Aromatic Diamines - Urine (UADM_I)
diamines = data %>%
  filter(names == "diamines") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,URX4TDA,URX6TDA,URX4MDA,URX5NDA,
         URXPPDA)
# URX4TDA: 2,4-Diaminotoluene (4TDA) (ng/mL)
# URX6TDA: 2,6-Diaminotoluene (6TDA) (ng/mL)
# URX4MDA: 4MDA (ng/mL)
# URX5NDA: 1,5-Diaminonaphthalene (5NDA) (ng/mL)
# URXPPDA: p-Phenylenediamine (PPDA) (ng/mL)


#### OUTCOMES ####
# Complete Blood Count with 5-Part Differential - Whole Blood (CBC_I)
blood_count = data %>%
  filter(names == "blood_count") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,LBXWBCSI,LBXLYPCT,LBXMOPCT,LBXNEPCT,
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


blood_pressure = data %>%
  filter(names == "blood_pressure") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,BPXSY1,BPXDI1,BPXSY2,BPXDI2,
         BPXSY3,BPXDI3,BPXSY4,BPXDI4)
# BPXSY1: Systolic: Blood pres (1st rdg) mm Hg
# BPXDI1: Diastolic: Blood pres (1st rdg) mm Hg
# BPXSY2: Systolic: Blood pres (2nd rdg) mm Hg
# BPXDI2: Diastolic: Blood pres (2nd rdg) mm Hg
# BPXSY3: Systolic: Blood pres (3rd rdg) mm Hg
# BPXDI3: Diastolic: Blood pres (3rd rdg) mm Hg
# BPXSY4: Systolic: Blood pres (4th rdg) mm Hg
# BPXDI4: Diastolic: Blood pres (4th rdg) mm Hg

body_measures = data %>%
  filter(names == "body_measures") %>%
  dplyr::select(data) %>%
  unnest() %>%
  #filter(SEQN %in% chem_nhanes$SEQN) %>%
  dplyr::select(SEQN,BMXWT,BMXRECUM,BMXHEAD,
         BMXHT,BMXBMI,BMXLEG,BMXARML,
         BMXARMC,BMXWAIST,BMDAVSAD)
# BMXWT: Weight (kg)
# BMXRECUM: Recumbent Length (cm)
# BMXHEAD: Head Circumference (cm)
# BMXHT: Standing Height (cm)
# BMXBMI: Body Mass Index (kg/m**2)
# BMXLEG: Upper Leg Length (cm)
# BMXARML: Upper Arm Length (cm)
# BMXARMC: Arm Circumference (cm)
# BMXWAIST: Waist Circumference (cm)
# BMDAVSAD: Average Sagittal Abdominal Diameter (cm)

# Possibility to add variables about DENTAL STATUS (either cathegorical or binary or teeth counts)
# https://wwwn.cdc.gov/nchs/nhanes/Search/DataPage.aspx?Component=Examination&CycleBeginYear=2015

#### JOIN #### 
df_chem = join_all(list(chromium_cobalt,copper_etc,cotinine,cotinine_urine, 
                   deet,diamines,ferritin,fluoride_plasma,fluoride_water,
                   folate,folate_serum,glycohemoglobin,iodine,mercury_blood,
                   mercury_urine,metals_blood,metals_urine,neonicotinoids,
                   pfas,phthalate,spec_arsernic_urine,steroid,trichomonas,
                   voc_blood,voc_urine), 
              by='SEQN', type='full')

df_out = join_all(list(blood_count,blood_pressure,
                       body_measures), 
                   by='SEQN', type='full')

df_cov = join_all(list(people,albumin_creatinine_nhanes,
                       apolipoprotein,bio,chol, 
                       arsenic,chlamydia), 
                   by='SEQN', type='full')

df = join_all(list(df_cov,
                   df_out,
                   df_chem), 
              by='SEQN', type='full')


#### SAVE ####
save(df, file = "~/SEM/data/nhanes_1516.RData")
save(df_chem, file = "~/SEM/data/nhanes_chem_1516.RData")
save(df_cov, file = "~/SEM/data/nhanes_cov_1516.RData")
save(df_out, file = "~/SEM/data/nhanes_out_1516.RData")


# Chemical specific data
df_metals = join_all(list(chromium_cobalt,copper_etc,mercury_blood,
                        mercury_urine,metals_blood,metals_urine,
                        spec_arsernic_urine), 
                   by='SEQN', type='full')
df_phalates_pfas = join_all(list(pfas,phthalate), 
                   by='SEQN', type='full')
save(df_metals, file = "~/SEM/data/nhanes_metals_1516.RData")
save(df_phalates_pfas, file = "~/SEM/data/nhanes_phalates_pfas_1516.RData")









