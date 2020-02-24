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

mysurvplot <- function(fit, group, data){
  fit %>%
    ggsurvplot_group_by(
      group.by = group, data = data,
      surv.median.line = "hv", conf.int = TRUE, pval = TRUE,
      risk.table = TRUE, risk.table.pos = "in", risk.table.fontsize = 3,
      ggtheme = theme_bw()
    ) %>%
    arrange_ggsurvplots() %>%
    print()
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
  select(kh, tte, status, Sex, Age, BMIc, WAISc, VM12c) %>%
  as.data.frame()

GGally::ggpairs(iki[, c("Age", "Sex", "tte")])

remove(list = ls()[!str_detect(ls(), "iki|my")])
haven::write_sas(iki, "iki.sas7bdat")

# non-parametric ----------------------------------------------------------

survfit(Surv(tte, status) ~ Age, data=iki) %>% mysurvplot("Sex", iki)
survfit(Surv(tte, status) ~ Sex, data=iki) %>% mysurvplot("Age", iki)

# semi-parametric ---------------------------------------------------------

## general
fit_cox <- coxph(Surv(tte, status) ~ Sex + Age + BMIc + WAISc + VM12c, data = iki) 
fit_cox %>% cox.zph() %>% ggcoxzph()
fit_cox %>% ggforest()

fit_cox %>% ggadjustedcurves("Sex"  , iki)
fit_cox %>% ggadjustedcurves("Age"  , iki)
fit_cox %>% ggadjustedcurves("WAISc", iki)
fit_cox %>% ggadjustedcurves("VM12c", iki)

# parametric --------------------------------------------------------------

survreg(Surv(tte, status) ~ Sex + Age + BMIc + WAISc + VM12c, data = iki) %>% plot()

# Imputation --------------------------------------------------------------

library(mice)

iki_miss <- iki_ok %>% 
  #  filter(tte > 0) %>%
  select(kh, tte, status, Sex, Age, bmi, WAIS, vm12) %T>%
  md.pattern(plot = TRUE, rotate.names = TRUE)
  
fit_mipo <- iki_miss %>% 
  mice(m=100, printFlag = FALSE) %>%
  with(coxph(Surv(tte, status) ~ Age + Sex + bmi + WAIS + vm12)) %>%
  pool(dfcom = nrow(iki_miss) - 6)


# forest plot -------------------------------------------------------------


library(forestplot)

fit_mipo %>% 
  summary() %>%
  mutate(p.value = format.pval(p.value, digits = 3, eps = 0.001),
         p.sig    = case_when(p.value < 0.001 ~ "***",
                              p.value < 0.01 ~ "**",
                              p.value < 0.05 ~ "*",
                              TRUE ~ ""),
         est   = exp(estimate),
         upper = exp(estimate + qnorm(0.975) * std.error),
         lower = exp(estimate + qnorm(0.025) * std.error)) %>%
  mutate(p.value  = paste0(p.value, p.sig),
         estimate = paste0(round(est, 3), "\n(", 
                           round(lower, 3), " - ", 
                           round(upper, 3), ")")) %>%
  select(term, estimate, p.value, est, upper, lower)-> x

forestplot(x[, 1:3], 
           txt_gp = fpTxtGp(label = list(gpar(cex = 0.7, fontfamily = "Times"))),
           title = "Hazard Ratio",
           mean  = x$est,
           lower = x$lower,
           upper = x$upper,
           zero   = 1,
           xlab   = "Relative Hazard Ratio",
           boxsize = 0.1,
           graph.pos  = 3)

ggforest(fit_cox)

library(broom); tidy(fit_cox); glance(fit_cox)

# References --------------------------------------------------------------
## https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
## https://stefvanbuuren.name/fimd/

## https://openaccess.leidenuniv.nl/bitstream/handle/1887/11456/01.pdf?sequence=6
## https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf

## https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html
## https://koppa.jyu.fi/avoimet/kirjasto/en/library-tutorial/citing-and-managing-references/citing-and-managing/how-to-cite