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
    tte    = elovuodet - Ika_tarkka,
    bmi    = paino / ((pituus/100)^2),
    BMI_   = u_BMIcat (paino, pituus),
    WAIS_  = u_WAIScat(WAIS),
    VM12_  = u_vmcat  (vm12, sukup, Ika),
    VM13c  = u_vmcat  (vm13, sukup, Ika),
    VM14c  = u_vmcat  (vm14, sukup, Ika),
    vm12   = vm12/kg2newton,
    vm13   = vm13/kg2newton,
    vm14   = vm14/kg2newton,
    Age    = factor(Ika, levels = c(75, 80)),
    Sex    = ifelse(sukup == "mies", "Male", "Female") %>% factor(levels =  c("Male", "Female")),
    dateEND= parse_date(substr(1000000 + if_else(is.na(kpvm), Haastpv, kpvm), 2, 9), format="%d%m%y"),
    dateSTX= parse_date(substr(1000000 + syntaika, 2, 9) , format="%d%m%y"), 
    statusX= ifelse(is.na(kpvm), 0 ,1)) %>% # mutate (tteX   = (dateEND - dateSTX)/365.25 + 100 - Age_tarkka)
  filter(!is.na(tte))

iki <- iki_ok %>% 
  select(kh, tte, status, Sex, Age, BMI_, WAIS_, VM12_) %>%
  as.data.frame()

GGally::ggpairs(iki[, c("Age", "Sex", "tte")])

remove(list = ls()[!str_detect(ls(), "iki|my")])

# foreign::write.foreign(iki, codefile = "iki",  datafile = "ikid", package = "SAS")

# non-parametric ----------------------------------------------------------

mysurvplot <- function(fit, group, data){
  fit %>%
    ggsurvplot_group_by(
      group.by = group, data = data,
      surv.median.line = "hv", conf.int = TRUE, 
      pval = TRUE, pval.size = 4, pval.coord = c(15, 0.8),
      risk.table = TRUE, risk.table.pos = "in", risk.table.fontsize = 3,
      ggtheme = theme_classic2(base_family = "Times"), palette = "Set2",
      font.family = "Times")
}
survfit(Surv(tte, status) ~ Age, data=iki) %>% mysurvplot("Sex", iki) -> p1
survfit(Surv(tte, status) ~ Sex, data=iki) %>% mysurvplot("Age", iki) -> p2
out_kmplots <- c(p1, p2)

# semi-parametric ---------------------------------------------------------

fit_cox  <- coxph(Surv(tte, status) ~ Sex + Age + BMI_ + WAIS_ + VM12_, data = iki)
fit_cox0 <- coxph(Surv(tte, status) ~ Sex + Age , data = iki) 
fit_cox %>% ggforest()

out_diag <- lapply(list(fit_cox0, fit_cox), . %>% cox.zph() %>% ggcoxzph(caption = "hi"))

out_adj <- lapply( list("Sex", "Age", "WAIS_", "VM12_"), 
                   function(x) ggadjustedcurves(
                     fit_cox, variable = x, data = iki,
                     ggtheme = theme_classic2(base_family = "Times"), palette = "Set2",
                     font.family = "Times") + labs(color = x) ) 

# Imputation --------------------------------------------------------------

library(mice)

if(!exists("n_mi")) n_mi <- 10

iki_miss <- iki_ok %>% 
  filter(tte > 0) %>%
  mutate(BMI_value = bmi, WAIS_value = WAIS, VM12_value = vm12) %>%
  select(kh, tte, status, Sex, Age, BMI_value, WAIS_value, VM12_value) %T>%
  md.pattern(plot = TRUE, rotate.names = TRUE)
  
fit_mipo <- iki_miss %>% 
  mice(m = n_mi, printFlag = FALSE) %>%
  with(coxph(Surv(tte, status) ~ Age + Sex + BMI_value + WAIS_value + VM12_value)) %>%
  pool(dfcom = nrow(iki_miss) - 6)


# forest plot -------------------------------------------------------------

library(ggforestplot)
library(broom)

m1 <- fit_cox0 %>% tidy() %>% mutate(Model = "Model 1")
m2 <- fit_cox %>% tidy() %>% mutate(Model = "Model 2")
m3 <- fit_mipo %>% summary() %>% mutate(Model = "Model 3") 
mm <- bind_rows(m1, m2, m3)

mm %>% 
  forestplot(name = term, pvalue = p.value, se = std.error,
             logodds = TRUE, colour = Model, shape = Model)


# estimation table --------------------------------------------------------

library(DT)
sketch = htmltools::withTags(table(
  class = 'display',
  thead(tr(th(rowspan = 2, 'Variable'),
           th(colspan = 2, 'Base Model'),
           th(colspan = 2, 'with Factorized Baseline Covariates'),
           th(colspan = 2, 'with Imputated Baseline Covariates')),
        tr(lapply(rep(c("Hazard Ratio (95% CI)", "p-value"), 3), th)))))

mm %>%
  mutate(len = ifelse((str_detect(term, "value")), 3, 3)) %>%
  transmute(term = term %>% 
              str_replace("BMI_value", "BMI (kg/m2)") %>%
              str_replace("WAIS_value", "WAIS Index") %>%
              str_replace("VM12_value", "Hand Grip Strength (kg)"),
            exp  = exp(estimate) %>% round(len),
            cil  = exp(estimate + qnorm(0.025) * std.error) %>% round(len),
            ciu  = exp(estimate + qnorm(0.975) * std.error) %>% round(len),
            pval = case_when(p.value<0.001 ~ "<0.001***", 
                             p.value<0.01  ~ paste0(round(p.value, 3), "**"),
                             p.value<0.05  ~ paste0(round(p.value, 3), "*"),
                             p.value<1.00  ~ paste0(round(p.value, 3))),
            ref  = as.factor(Model)) %>%
  transmute(term = str_replace(term, "BMI_|WAIS_|VM12_", "-- "),
            est  = paste0(exp, "<br/>(", cil, "-", ciu, ")"),
            pval = pval,
            ref = ref) -> tmp


filter(tmp, ref == levels(ref)[1]) %>% 
  full_join(filter(tmp, ref == levels(ref)[2]), by = "term") %>%
  full_join(filter(tmp, ref == levels(ref)[3]), by = "term") %>%
  add_row(term = "SexMale", est.x = "Ref.", est.y = "Ref.", est = "Ref.", .before = 1) %>%
  add_row(term = "Age75",   est.x = "Ref.", est.y = "Ref.", est = "Ref.", .before = 3) %>%
  add_row(term = "BMI_Missing", est.y = "Ref.", .before = 5) %>%
  add_row(term = "WAIS_Missing", est.y = "Ref.", .before = 10) %>%
  add_row(term = "VM12_Missing", est.y = "Ref.", .before = 16) %>%
  add_row(term = "Continuous Variable", .before = 22) %>%
  select(-ref.x, -ref.y, -ref) %>%
  datatable(rownames = FALSE, container = sketch, escape = FALSE,
            options = list(dom="t", pageLength=nrow(.))) %>%
  formatStyle(2:10, textAlign = 'center') %>%
  formatStyle(1:10, target = 'row', 
              backgroundColor = styleEqual(c("Ref.", "Continuous Variable"), rep("lightgrey",2))) %>%
  formatStyle(c(1,3,5), `border-right` = "solid 1px") -> out_est_dt




# References --------------------------------------------------------------
## https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
## https://stefvanbuuren.name/fimd/

## https://openaccess.leidenuniv.nl/bitstream/handle/1887/11456/01.pdf?sequence=6
## https://rpkgs.datanovia.com/survminer/survminer_cheatsheet.pdf

## https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html
## https://koppa.jyu.fi/avoimet/kirjasto/en/library-tutorial/citing-and-managing-references/citing-and-managing/how-to-cite
