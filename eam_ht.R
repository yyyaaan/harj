library(survival)
library(survminer)
library(tidyverse)

load("eam_0.RData"); remove(list = ls()[!str_detect(ls(), "iki")])



# Supporting utilities ----------------------------------------------------

kg2newton <<- 9.80665
factorcls <<- c("Missing", paste("Class", 1:5))
factorbmi <<- c("Missing", "Underweight", "Normal weight", "Overweight", "Obese")
factorAge <<- c(75, 80)
factorsuk <<- c("mies", "nainen")

u_vmcat <- function(vm, Sex, Age){
  if(is.na(vm)) return(factor("Missing", factorcls))
  vm  <- vm/kg2newton
  out <- factor("Class 1", factorcls)
  if(Sex == "nainen" && Age == 75){
    if(vm >= 15.8) out <- factor("Class 2", factorcls)
    if(vm >= 19.7) out <- factor("Class 3", factorcls)
    if(vm >= 22.7) out <- factor("Class 4", factorcls)
    if(vm >= 26.1) out <- factor("Class 5", factorcls)
  }
  if(Sex == "nainen" && Age == 80){
    if(vm >= 11.3) out <- factor("Class 2", factorcls)
    if(vm >= 14.7) out <- factor("Class 3", factorcls)
    if(vm >= 17.5) out <- factor("Class 4", factorcls)
    if(vm >= 21.1) out <- factor("Class 5", factorcls)
  }
  if(Sex == "mies" && Age == 75){
    if(vm >= 28.6) out <- factor("Class 2", factorcls)
    if(vm >= 33.5) out <- factor("Class 3", factorcls)
    if(vm >= 39.0) out <- factor("Class 4", factorcls)
    if(vm >= 43.4) out <- factor("Class 5", factorcls)
  }
  if(Sex == "mies" && Age == 80){
    if(vm >= 20.3) out <- factor("Class 2", factorcls)
    if(vm >= 24.7) out <- factor("Class 3", factorcls)
    if(vm >= 31.0) out <- factor("Class 4", factorcls)
    if(vm >= 34.9) out <- factor("Class 5", factorcls)
  }
  return(out)
}

u_BMIcat <- function(weight, height){
  if(is.na(weight) || is.na(height)) return(factor("Missing", factorbmi))
  h <- ifelse(height > 3, height / 100, height)
  bmi <- weight / (h^2)
  out <- factor("Underweight", factorbmi)
  if (bmi >= 18.5) out <- factor("Normal weight", factorbmi)
  if (bmi >= 25.0) out <- factor("Overweight", factorbmi)
  if (bmi >  30.0) out <- factor("Obese", factorbmi)
  return(out)
}

u_WAIScat <- function(wais){
  if(is.na(wais)) return(factor("Missing", factorcls))
  out <- factor("Class 1", factorcls)
  if(wais >= 12) out <- factor("Class 2", factorcls)
  if(wais >= 17) out <- factor("Class 3", factorcls)
  if(wais >= 21) out <- factor("Class 4", factorcls)
  if(wais >= 28) out <- factor("Class 5", factorcls)
  return(out)
}

# data --------------------------------------------------------------------

iki_ok <- ikivihreat %>% 
  rowwise() %>% 
  mutate( 
    status = ifelse(Elostatus == "Kuollut", 1, 0),
    tte    = elovuodet - Ika_tarkka ,
    bmi    = paino / ((pituus/100)^2),
    BMIc   = u_BMIcat (paino, pituus),
    WAISc  = u_WAIScat(WAIS),
    VM12c  = u_vmcat  (vm12, sukup, Ika),
    VM13c  = u_vmcat  (vm13, sukup, Ika),
    VM14c  = u_vmcat  (vm14, sukup, Ika),
    Age    = factor(Ika, levels = c(75, 80)),
    Sex    = ifelse(sukup == "mies", "Male", "Female") %>% factor(levels =  c("Male", "Female")),
    dateEND= parse_date(substr(1000000 + if_else(is.na(kpvm), Haastpv, kpvm), 2, 9), format="%d%m%y"),
    dateSTX= parse_date(substr(1000000 + syntaika, 2, 9) , format="%d%m%y"), 
    statusX= ifelse(is.na(kpvm), 0 ,1)) # %>% mutate (tteX   = (dateEND - dateSTX)/365.25 + 100 - Age_tarkka)

iki <- iki_ok %>% 
  filter(tte > 0) %>%
  select(kh, tte, status, Sex, Age, BMIc, WAISc, VM12c, VM13c, VM14c) %>%
  as.data.frame()

GGally::ggpairs(iki[, c("Age", "Sex", "tte")])

## clean up environment
remove(list = ls()[!str_detect(ls(), "iki|my")])



# non-parametric ----------------------------------------------------------

survfit(Surv(tte, status) ~ Age, data=iki) %>% mysurvplot("Sex", iki)
survfit(Surv(tte, status) ~ Sex, data=iki) %>% mysurvplot("Age", iki)

# semi-parametric ---------------------------------------------------------

## general
fit_cox <- coxph(Surv(tte, status) ~ Sex + Age + BMIc + WAISc + VM12c, data = iki) 
fit_cox %>% cox.zph() %>% ggcoxzph()
fit_cox %>% ggforest()

fit_cox %>% ggadjustedcurves(data = iki, variable = "Sex") 
fit_cox %>% ggadjustedcurves(data = iki, variable = "Age")
fit_cox %>% ggadjustedcurves(data = iki, variable = "BMIc")
fit_cox %>% ggadjustedcurves(data = iki, variable = "WAISc")
fit_cox %>% ggadjustedcurves(data = iki, variable = "VM12c")


## stratified
for(stratum in levels(iki$Age)){
  ikis <- iki %>% filter(Age == stratum)
  coxph(Surv(tte, status) ~ Sex + BMIc + WAISc + VM12c, data = ikis) %>% 
    ggforest(main = paste("Hazard Ratio | Stratum: Age =", stratum), fontsize = 0.9) %>%
    print()
}

# parametric --------------------------------------------------------------

survreg(Surv(tte, status) ~ Sex + Age + BMIc + WAISc + VM12c, data = iki) %>% plot()
  ggsurvplot()

# Imputation --------------------------------------------------------------

library(mice)
iki_miss <- iki_ok %>% 
  select(kh, tte, status, Sex, Age, bmi, WAIS, vm12) %T>%
  md.pattern(plot = TRUE, rotate.names = TRUE)
  
iki_miss %>% 
  mice(m=1000, printFlag = FALSE) %>%
  with(coxph(Surv(tte, status) ~ Age + Sex + bmi + WAIS + vm12)) %>%
  pool(dfcom = nrow(iki_miss) - 6) -> fit_mipo
  
## stratified
for(stratum in levels(iki$Age)){
  
  cat("\n=========\nAge Group ", stratum, "\n=========\n")
  ikis <- iki_miss %>% filter(Age == stratum) 
  ikis %>% 
    mice(m = 100, printFlag = FALSE) %>%
    with(coxph(Surv(tte, status) ~ Sex + bmi + WAIS + vm12)) %>%
    pool(dfcom = nrow(ikis) - 5) %>%
    summary(conf.int = TRUE) %>% print()
}



# playground --------------------------------------------------------------

library(forestplot)

fit_cox %>% 
  summary() %>% 
  extract2(7) %>%
  data.frame() %>%
  transmute(rowname = rownames(.) ,
            value   = rownames(.) %>% str_remove("Sex|Age|BMIc|WAISc|VM12c"),
            Estimate = paste0(round(exp(coef), 3), "\n(", 
                              round(exp(coef + qnorm(0.025) * se.coef.),3),
                              "-", round(exp(coef + qnorm(0.975) * se.coef.),3),")"),
            p.value  = ifelse(Pr...z.. <0.001, "<.001", round(Pr...z..,  3)),
            p.sig    = case_when(Pr...z.. < 0.001 ~ "***",
                                 Pr...z.. < 0.01 ~ "**",
                                 Pr...z.. < 0.05 ~ "*",
                                 TRUE ~ ""),
            est   = exp(coef),
            upper = exp(coef + qnorm(0.975) * se.coef.),
            lower = exp(coef + qnorm(0.025) * se.coef.)) %>%
  mutate(p.value = paste0(p.value, p.sig),
         variable= str_remove(rowname, value)) -> x; x

iki %>% 
  reshape2::melt(id.vars = c("kh", "tte", "status")) %>%
  group_by(variable, value) %>%
  summarise(n()) %>%
  left_join(x) -> aaa
  
x$rowname 
# Rowname = row.names(.) %>%
#   str_replace("SexFemale", "Sex Female vs. Male") %>%
#   str_replace("BMIc", "BMIc ") %>%
#   str_replace("cClass", "c Class") %>%
#   ifelse(str_detect(., ":"), paste(., " vs. Missing"), .)  %>% 
#   str_replace("Age80", "Age 80 vs. 75"),


gpl <- gpar(lty=1)
cls <- fpColors(box="royalblue",line="darkblue", summary="royalblue")
txt <- fpTxtGp(label = list(gpar(cex = 0.7, fontfamily = "Times")))
forestplot(x[, 1:3], txt_gp = txt,
           mean  = x$est,
           lower = x$lower,
           upper = x$upper,
           zero   = 1,
           clip   = c(0, 3.3),
           xlab   = "Relative Hazard Ratio",
           boxsize = 0.3,
           graph.pos  = 3,
           col = cls,
           hrzl_lines = list("2"=gpl, "3"=gpl, "7"=gpl, "12"=gpl))


ggforest(fit_cox)

# References --------------------------------------------------------------
## https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
## https://stefvanbuuren.name/fimd/

## https://openaccess.leidenuniv.nl/bitstream/handle/1887/11456/01.pdf?sequence=6
## https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf
