
# clr 101 -----------------------------------------------------------------

options(stringsAsFactors = T)
aineisto=read.table("altruismi.dat")
str(aineisto)
aineisto$B18E[aineisto$B18E == 133.9] <- 13.9

require(cluster)
dis <- daisy(aineisto)
summary(dis)

# ylm 608 -----------------------------------------------------------------
new <- data.frame(vuosi = c(2015:2020, 2015:2020),
                  sukup = c(rep(0,6), rep(1,6)),
                  vaki = 1e5)
fitted <- predict(fit607, newdata = new, type = "response", se.fit = TRUE)
data.frame(vuosi = new$vuosi,
           sukup = new$sukup,
           est = fitted$fit,
           lwr = fitted$fit - 1.96 * fitted$se.fit,
           upr = fitted$fit + 1.96 * fitted$se.fit)


# ylm 607 -----------------------------------------------------------------

vk <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/verenkierto.txt",header=TRUE)

fit607 <- glm(kuolleet ~ vuosi * sukup + offset(log(vaki)), 
              family = poisson, data = vk)
cbind(exp(coef(fit607)), exp(confint(fit607)))

b <- coef(fit607)
vk$lk <- vk$kuolleet / vk$vaki * 100000
plot(lk ~ vuosi, data = vk[vk$sukup == 1,], col = "red", ylim = c(25,350),
     main = "circulation disease death rate, with fit")
points(lk ~vuosi, data = vk[vk$sukup == 0,], col = "blue")
curve(exp(b[1] + b[2] * x + b[3] + b[4] * x) * 1e5, col = "red",
      from = 1970, to = 2014, add = TRUE)
curve(exp(b[1] + b[2] * x ) * 1e5, col = "blue",
      from = 1970, to = 2014, add = TRUE)
legend("topright", c("man","woman"), col = c("blue","red"), pch = 1)

# ylm 606 -----------------------------------------------------------------

library(MASS)
fit606 <- glm.nb(Colias ~ time + habitat + building, data = perhoset)
cbind(exp(coef(fit606)), exp(confint(fit606)))


# ylm 605 -----------------------------------------------------------------
  #ref
F_testi <- function(K, X, beta, m, sigma)
{
  q <- qr(K)$rank; n <- nrow(X); p <- ncol(X)
  C <- solve(t(K) %*% solve(t(X) %*% X) %*% K)
  F <- t(t(K) %*% beta - m) %*% C %*% (t(K) %*% beta - m) / (q * sigma)
  p_arvo <- 1 - pf(F, q, n - p - 1)
  list(F = F, p_arvo = p_arvo)
}

K <- diag(nrow = 3)
m <- matrix(rep(0,3), nrow = 3)
X <- data.frame(x1 = as.numeric(perhoset$habitat == "Mixed"),
                x2 = as.numeric(perhoset$habitat == "Short"),
                x3 = as.numeric(perhoset$habitat == "Tall"))
X <- as.matrix(X)
sigma2 <- t(residuals(fit603)) %*% residuals(fit603) / (n - length(b) -1)
beta <- matrix(b[3:5], nrow = 3)

F_testi(K, X, beta, m, sqrt(sigma2))


# ylm 603 -----------------------------------------------------------------

options(stringsAsFactors = TRUE)
perhoset <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/perhoset.txt",header=TRUE)

  # no site in model
fit603 <- glm(Colias ~ time + habitat + building, 
              family = poisson, 
              data = perhoset)
cbind(exp(coef(fit603)), exp(confint(fit603)))

b <- coef(fit603)
par(mfrow = c(2,2))
  #marginal plot on time
plot(Colias ~ time, pch = 18, col = "grey", 
     data = perhoset[perhoset$habitat == "Hayfield",], 
     main = "Yearly Change: Hayfield Habitat")
curve(exp(b[1] + b[2] * x + b[6] * mean(perhoset$building)), 
      from = 1, to = 5, add = TRUE)
plot(Colias ~ time, pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Mixed",], 
     main = "Yearly Change: Mixed Habitat")
curve(exp(b[1] + b[2] * x  + b[3] + b[6] * mean(perhoset$building)), 
      from = 1, to = 5, add = TRUE)
plot(Colias ~ time,  pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Short",], 
     main = "Yearly Change: Short Habitat")
curve(exp(b[1] + b[2] * x + b[4] + b[6] * mean(perhoset$building)), 
      from = 1, to = 5, add = TRUE)
plot(Colias ~ time,  pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Tall",], 
     main = "Yearly Change: Tall Habitat")
curve(exp(b[1] + b[2] * x + b[5] + b[6] * mean(perhoset$building)), 
      from = 1, to = 5, add = TRUE)

  # marginal plot on building
plot(Colias ~ building, pch = 18, col = "grey", 
     data = perhoset[perhoset$habitat == "Hayfield",], 
     main = "Effect of Building Density: Hayfield Habitat")
curve(exp(b[1] + b[2] * 3 + b[6] * x), 
      from = 0, to = 20, add = TRUE)
plot(Colias ~ building, pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Mixed",], 
     main = "Effect of Building Density: Mixed Habitat")
curve(exp(b[1] + b[2] * 3  + b[3] + b[6] * x), 
      from = 0, to = 20, add = TRUE)
plot(Colias ~ building,  pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Short",], 
     main = "Effect of Building Density: Short Habitat")
curve(exp(b[1] + b[2] * 3 + b[4] + b[6] * x), 
      from = 0, to = 20, add = TRUE)
plot(Colias ~ building,  pch = 18, col = "grey",
     data = perhoset[perhoset$habitat == "Tall",], 
     main = "Effect of Building Density: Tall Habitat")
curve(exp(b[1] + b[2] * 3 + b[5] + b[6] * x), 
      from = 0, to = 20, add = TRUE)


# ylm 604 -----------------------------------------------------------------




# ylm 505 -----------------------------------------------------------------

lapset <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/lapset.txt")
  #define new vars
lapset$lapsia <- 0
lapset$lapsia[which(lapset$parit == 1)] <- 1
lapset$lapsia[which(lapset$parit > 1)] <- 2
lapset$kesk <- ifelse(lapset$syntpaino < 2500, 1, 0)
head(lapset)


fit1 <- glm(kesk ~ as.factor(tup) * as.factor(lapsia) + as.factor(sukup) ,
            family = binomial, data = lapset)
beta <- cbind(est = coef(fit1), confint(fit1))
or <- exp(beta)
colnames(or) <- c("odds.ratio", "OR 2.5%", "OR 97.5%")
or

beta[2,] + beta[4,] + beta[7,]

fit2 <- glm(kesk ~ as.factor(tup) * as.factor(lapsia) + as.factor(sukup) + rviik,
            family = binomial, data = lapset)
exp(cbind(odds.ratio = coef(fit2), confint(fit2)))

# ylm 504 -----------------------------------------------------------------

tul <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/korvatul.txt",header=T)
# vaste ja design
y <- tul[,1]
X0 <- rep(1,length(y))
X1 <- tul[,2]
X2 <- 1 * (tul[,3] == 2)
X3 <- 1 * (tul[,3] == 3)
X <- cbind(X0, X1, X2, X3)
# Alkuarvo betalle ja epsilonille
beta <- c(0.1, 0.1, 0.1, 0.1)
eps <- 
  
  while (eps > 1e-4) {
    # eta, mu ja score-funktio
    eta <-
    pi <- exp(eta) / (1 + exp(eta))
    S <-
    # W matrix
    Wii <-
    W <- diag(c(Wii))
    # Fisherin informaatiomatriisi
    I <-
    # Iteraatioaskel
    betanew <- beta + solve(I) %*% S
    # eps pysäytyskriteeriä varten
    eps <- t(betanew-beta) %*% (betanew-beta)
    # päivitetään \beta
    beta <- betanew
  }
# Listataan SU-estimaatit ja keskivirheet
list(beta = , se = )


# put 504 -----------------------------------------------------------------

# pars = c( mu1, mu2, sigma1, sigma2, lambda)
likely <- function(x, pars){
  n <- nrow(x)
  p <- ncol(x)
  mu1 <- matrix(pars[1:p], ncol = p)
  mu2 <- matrix(pars[(p+1):(2*p)], ncol = p)
  sigma1 <- matrix(pars[(2*p + 1): (2*p + p*p)], ncol = p)
  sigma2 <- matrix(pars[(2*p + p*p + 1): (2*p + 2* p*p)], ncol = p)
  lambda <- pars[-1]
  
  z <- rbinom(n, 1, lambda) # simulate z \in {0, 1}
  logf <- 0
  for (i in 1:n) {
    logf <- logf + 0.5 * p * log(2 * pi)
      ifelse(z == 1, 
             log(lambda) - 0.5 * log(det(sigma1)) - 
               0.5 * (x[i,] - mu1) %*% solve(sigma1) %*% t(x[i,] - mu1),
             log(1 - lambda) - 0.5 * log(det(sigma2)) - 
               0.5 * (x[i,] - mu2) %*% solve(sigma2) %*% t(x[i,] - mu2))
      
  }
  logf
}


# library(MASS)
# xa <- mvrnorm(n = 1e3, mu = c(0,1), Sigma = matrix(c(1,0.2,0.2,2), ncol = 2))
# xb <- mvrnorm(n = 1e3, mu = c(8,9), Sigma = matrix(c(5,0.2,0.2,9), ncol = 2))
# x <- xa + xb 
# 
# init <- c(mu1 = matrix(c(1,2), ncol = 2),
#           mu2 = matrix(c(8,9), ncol = 2),
#           sigma1 = matrix(c(1,0.5,0.5,1), ncol = 2),
#           sigma2 = matrix(c(3,0.5,0.5,3), ncol = 2),
#           lambda = 0.3) # will always be numeric vector
# 
# optim(init, fn = function(pars) likely(x,pars))



# put 503 -----------------------------------------------------------------

#folllow chapter 8.5 normaalijakauman sekoitus

library(MASS)

xa <- mvrnorm(n = 1e3, mu = c(-5,-5), Sigma = matrix(c(1,0,0,1), ncol = 2))
xb <- mvrnorm(n = 1e3, mu = c(5,5), Sigma = matrix(c(1,0,0,1), ncol = 2))
p <- rbinom(1e3, 1, 0.7)
x <- xa
x[which(p==0),] <- xb [which(p==0),]
plot(x)

  #init value 
init <- list(mu1 = matrix(c(1,2), ncol = 2),
             mu2 = matrix(c(8,9), ncol = 2),
             sigma1 = matrix(c(1,2,3.9,1), ncol = 2),
             sigma2 = matrix(c(3,0.5,0.9,3), ncol = 2),
             lambda = 0.3)

myEM <- function(init, x, maxIte) {
  p <- length(init$mu1) # dimension of normal
  theta <- list()
  theta[[1]] <- init
  
  for (i in 2:maxIte) {
    print(i)
    theta.p<- theta[[i-1]]
    
    R1 <- numeric(nrow(x))
    R2 <- numeric(nrow(x))
    # use log for accuracy
    for (j in 1:nrow(x)) {
      logf1 <-  -0.5 * log( det(theta.p$sigma1)) - 
        0.5 * (x[j,] - theta.p$mu1) %*% solve(theta.p$sigma1) %*% t(x[j,] - theta.p$mu1) +
        p / 2 * log(2 * pi)
      logf2 <-  -0.5 * log( det(theta.p$sigma2)) - 
        0.5 * (x[j,] - theta.p$mu2) %*% solve(theta.p$sigma2) %*% t(x[j,] - theta.p$mu2) +
        p / 2 * log(2 * pi)
      logRa <- log(theta.p$lambda) + logf1
      logRb <- log(1 - theta.p$lambda) + logf2
      R1[j] <- exp( logRa - log(exp(logRa) + exp(logRb)))
      R2[j] <- 1- R1[j]
    }
    
    lambda <- sum(R1) / length(R1)
      # point-wise multiplication
    mu1 <- colSums(R1 * x) / sum(R1)
    mu2 <- colSums(R2 * x) / sum(R2)
      # calculate Rt(xi - mu1)(xi - mu1) !new mu
    rxx1 <- 0; rxx2 <-0
    for (j in 1: nrow(x)) {
      rxx1 <- rxx1 + R1[i] * t(x[2,] - mu1) %*% (x[2,] - mu1)
      rxx2 <- rxx2 + R2[i] * t(x[2,] - mu2) %*% (x[2,] - mu2)
    }
    
    theta[[i]] <- list(mu1 = mu1,
                       mu2 = mu2,
                       sigma1 = rxx1 / sum(R1),
                       sigma2 = rxx2 / sum(R2),
                       lambda = lambda)
  }
  
  return(theta)
}


# ylm 406 -----------------------------------------------------------------

korvatul <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/korvatul.txt")
fit <- glm(tulehdus ~ tupak * factor(lapsia), family = binomial, data = korvatul)
summary(fit)

table(korvatul$lapsia) / nrow(korvatul)

# ylm 402 -----------------------------------------------------------------

Sc <- function(y,mu,n) {
  n * (mean(y) - mu) 
}
Sc.p <- function(y, lambda){
  - length(y) + sum(y) / lambda
}
Sc.b <- function(y, pi){
  sum(y) / pi - (length(y) - sum(y)) / (1-pi)
}

n <- 10
par(mfrow= c(1,3) )

y <- rnorm(n, 4, 1); mu <- seq(0,6,0.1)
plot(mu, Sc(y, mu, n), type = "l", xlab = "mu", ylab = "S(mu)", main = "N(4,1)")
abline(0,0)
for (i in 1:19) {
  y <- rnorm(n, 4, 1)
  lines(mu, Sc(y, mu, n), type = "l")
}

y <- rpois(n, 4); lambda <- seq(1, 10, 0.1)
plot(lambda, Sc.p(y, lambda), type = "l", xlab = "lambda", ylab = "S(lambda)", main = "Poisson(4)")
abline(0,0)
for (i in 1:19) {
  y <- rpois(n, 4)
  lines(lambda, Sc.p(y, lambda), type = "l")
}

y <- rbinom(n, 1, 0.4); pi <- seq(0,1,0.01)
plot(pi, Sc.b(y, pi), type = "l", xlab = "pi", ylab = "S(pi)", main = "Bin(1,0.4)")
abline(0,0)
for (i in 1:19) {
  y <- rbinom(n, 1, 0.4)
  lines(pi, Sc.b(y, pi), type = "l")
}

# ylm 307 -----------------------------------------------------------------

# includes the interactive items
fith <- lme(sciescore ~ matem * motiv + gender + ESCS + urban + region,
            random = ~ 1 | SCHOOLID)
summary(fith)
anova.lme(fit2, fith)
plot(random.effects(fith))



# ylm 306 -----------------------------------------------------------------
  #a) read data
pisa <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/pisafull.txt")
  #b) ggpairs plot; this gives the most useful plots for explore data
library(GGally)
ggpairs(pisa)
  #c) d)
attach(pisa)
library(nlme)
  #e) generalized linear model fitting
fit <- gls(sciescore ~ matem + gender + ESCS + motiv + urban + region,
           correlation = corCompSymm(form = ~ 1 | SCHOOLID))
summary(fit)
  #f) 
s <- vcov(fit)[7:9, 7:9]
beta <- matrix(coef(fit)[7:9], nrow = 3)
F <- t(beta) %*% solve(s) %*% beta / 3
qF <- qf(.95, df1 = 3, df2 = nrow(pisa) - length(coef(fit)) - 1)
cat("F is caculated to be", F, "\nF-distribution quantile is", qF)
  #g)
fit2 <- gls(sciescore ~ matem * motiv + gender + ESCS + urban + region,
            correlation = corCompSymm(form = ~ 1 | SCHOOLID))
summary(fit2)
  #h)
df <- data.frame(school = SCHOOLID, residual = residuals(fit2))
ggplot(df, aes(school, residual)) + geom_jitter()
  #i) glm is indeed better
fit.lm <- lm(sciescore ~ matem * motiv + gender + ESCS + urban + region)
anova.lme(fit.lm, fit2)
  #j)
fit2.ml <- gls(sciescore ~ matem * motiv + gender + ESCS + urban + region,
            correlation = corCompSymm(form = ~ 1 | SCHOOLID), method = "ML")
cbind("ML - se" = coef(summary(fit2.ml))[,2],
      "REML - se" = coef(summary(fit2))[,2])

# ylm 304 -----------------------------------------------------------------

pisa <- read.table("http://users.jyu.fi/~slahola/files/glm2_datoja/pisa.txt", header = T)
mat <- pisa$matem
sp <- as.numeric(pisa$sukup == "tytto")
sij <- as.numeric(pisa$koulusij == "maaseutu")
ita <- as.numeric(pisa$koulualue == "Ita-Suomi")
lansi <- as.numeric(pisa$koulualue == "Lansi-Suomi")
pohjoinen <- as.numeric(pisa$koulualue == "Pohjois-Suomi")
X <- cbind(rep(1,200), mat, sp, sij, ita, lansi, pohjoinen)
y <- pisa$mpist

# Lasketaan SU-estimaatti

betahat <- solve(t(X) %*% X) %*% t(X) %*% y
betahat

# Mallin keskivirheen laskemiseksi tarvitaan jäännökset

resid <- y - X %*% betahat

# Mallin keskivirhe

n <- length(y)
p <- ncol(X) - 1
sigmahat <- (sqrt(t(resid) %*% resid / (n - p - 1)))[1,1]
sigmahat

# Lasketaan sitten regressiokerrointen keskivirheet

covbeta <- sigmahat^2 * solve(t(X)%*%X)
sdbeta <- sqrt(diag(covbeta))
round(sdbeta,2)

# Regressiokerrointen 95% luottamusvälit

ala <- betahat - qt(0.975, n - p - 1) * sdbeta
yla <- betahat + qt(0.975, n - p - 1) * sdbeta
cbind(ala,yla)

# TEHTÄVÄ: mitä luottamusvälit kertovat sinulle regressiokerrointen 
# merkitsevyydestä?

# TEHTÄVÄ: millä R:n funktiolla saat luottamusvälit laskettua 
# automaattisesti? 

# Selitysasteen laskemiseksi tarvitsemme kokonaisvaihtelun

ym <- y - mean(y)
SST <- t(ym) %*% ym 

# ja jäännösvaihtelun

SSE <- t(resid) %*% resid

# Selitysaste on silloin 

1 - SSE/SST

# TEHTÄVÄ: mitä selitysaste kertoo?

# Testataan sitten sukupuoleen liittyvää regressiokerrointa, 
# H_0: \beta_2 = 0

# t-testisuure

T <- betahat[3]/sdbeta[3]
T

# Kaksisuuntaisen testin kriittinen arvo on nyt

qt(0.025, n - p - 1)

# Nollahypoteesiä ei siis voida hylätä merkitsevyystasolla 0.05. 
# Voimme toki laskea tarkan p-arvonkin

2 * (1 - pt(abs(T), n - p - 1))

# Kirjoitetaan funktio, joka laskee F-testisuureen arvon
# yleiselle lineaariselle hypoteesille

F_testi <- function(K, X, beta, m, sigma)
{
  q <- qr(K)$rank; n <- nrow(X); p <- ncol(X)
  C <-  solve(t(K) %*% solve(t(X) %*% X) %*% K)
  F <- t(t(K) %*% beta - m) %*% C %*% (t(K) %*% beta - m) / (q * sigma^2)
  p_arvo <- 1 - pf(F, q, n - p - 1)
  list(F = F, p_arvo = p_arvo)
}

# TEHTÄVÄ: lue luentomonisteesta, mitä eroa on F-testillä ja Waldin testillä, 
# joka näyttää lähes samalta..

# Testataan F-testillä, tarvitaanko yhtään selittäjää, 
# H_0: \beta_1 = ... = \beta_6 = 0

K <- diag(7); K <- K[,-1]
K
F_testi(K, X, beta = betahat, m = rep(0,6), sigma = sigmahat)

# TEHTÄVÄ: Tarkista, että K^T beta = m todella johtaa haluamaamme hypoteesiin.

# Nollahypoteesi voidaan hylätä merkitsevyystasolla 0.05.

# Testataan, tarvitaanko koulualuetta malliin mukaan, 
# H_0: \beta_4 = \beta_5 = \beta_6 = 0

K <- matrix(rep(0,7*3),7,3); K[5,1] <- K[6,2] <- K[7,3] <- 1
K
F_testi(K, X, beta = betahat, m = rep(0,3), sigma = sigmahat)

# Nollahypoteesi voidaan hylätä merkitsevyystasolla 0.05.

# Testataan, ovatko Itä- ja Pohjoissuomeen liittyvät 
# regressiokertoimet yhtäsuuret, # H_0: \beta_4 = \beta_6 

K <- matrix(c(0,0,0,0,1,0,-1),7,1)
K
F_testi(K, X, beta = betahat, m = 0, sigma = sigmahat)

# Nollahypoteesiä ei voida hylätä merkitsevyystasolla 0.05.

# Lasketaan R:n predict-funktion avulla sovitteet 
# esimerkiksi Etelä-Suomen kaupunkikoulujen pojille,
# joilla matematiikan arvosana on 5-10:

newdatap <- data.frame(matem = seq(5,10), sukup = "poika",
                       koulusij = "kaupunki", koulualue = "Etela-Suomi")
newdatap

fit <- lm(mpist ~ matem + sukup + koulusij + koulualue,
          data = pisa)
predict(fit, newdatap)

# TEHTÄVÄ: ennusta vastaavat pisteet Länsi-Suomen tytöille.

# Tehdään sitten inferenssiä kahdella erilaisella bootstrap-menetelmällä

# Lasketaan ensin luottamusvälit parametrittomalla bootstrap-menetelmällä

bootstrap <- function(y,X) 
{
  n <- length(y)
  s <- sample(1:n, replace=TRUE)
  out <- lm(y[s] ~ X[s,] - 1)
  coef(out)
}

# TEHTÄVÄ: Käytä ensin yhtä bootstrap-otosta M <- 1 ja katso, mitä b.boot sisältää.
# Lisää sitten otosten määrää M <- 10000


b.boot <- replicate(M, bootstrap(y,X))
tab.boot <- apply(b.boot, 1, quantile, c(.025, .975))
colnames(tab.boot) <- colnames(X)
t(tab.boot)

# TEHTÄVÄ: varmista, että ymmärrät, millä yllä oleva koodi laskee.

# TEHTÄVÄ: tutki, millainen vaikutus bootstrap-otosten määrällä on tuloksiin. 
# a) M <- 10, b) M <- 100, c) M <- 1000

# Jos mallioletusten voidaan olettaa olevan voimassa, voidaan käyttää  
# parametrista bootstrap-menetelmää. Tämän jätämme ensi viikon 
# harjoitustehtäväksi

# Testataan vielä sukupuoleen liittyvää regressiokerrointa
# parametrittoman bootstrap-testin avulla

bootstrap <- function(y,X) 
{
  n <- length(y)
  s <- sample(1:n, replace=TRUE)
  out <- lm(y[s] ~ X[s,] - 1)
  coef(summary(out))[3,1:2]
}

# TEHTÄVÄ: Käytä ensin yhtä bootstrap-otosta M <- 1 ja katso, mitä b.boot sisältää.
# Lisää sitten otosten määrää M <- 10000

M <- 1
b.boot <- replicate(M, bootstrap(y,X))
z <- (b.boot[1,]-betahat[3])/b.boot[2,] 

# TEHTÄVÄ: plottaa z-arvojen jakauma 
plot(density(z)) 

# p-arvo kaksisuuntaiselle vastahypoteesille
1 - 2 * sum(abs(betahat[3] / sdbeta[3]) > z) / M

# TEHTÄVÄ: varmista, että ymmärrät, millä yllä oleva koodi laskee.

# TEHTÄVÄ: bootstrappaaminen ei vaadi taoreettisia jakaumatuloksia. Pohdi, miksi ei siis
# aina käytettäisi bootstrappausta?


# sda 202 -----------------------------------------------------------------

fun <- function(r) exp(-5 * r) - 1
s <- matrix(c(0.2, 0.5, 0.4, 0.9, 0.2, 0.1, 0.5, 0.8), ncol=2)
y <- matrix(c(12,6,10,2),ncol=1)
s0 <- matrix(c(0.5,0.4),ncol=1)
u <- y - 5.0
C <- fun(as.matrix(dist(s,diag=TRUE)))
s.aug <- matrix(c(s[,1], s0[1,1], s[,2], s0[2,1]),ncol=2)
cs <- as.matrix( fun(as.matrix(dist(s.aug))[5,1:4]))
lambda0 <- solve(C) %*% cs
hat.u <- t(lambda0) %*% u
hat.y <- hat.u + 5 
k.var <- 1 - t(cs) %*% solve(C) %*% cs


library(geoR)
y.geo <- as.geodata(matrix(c(s[,1], s[,2], y[,1]), ncol = 3))
simplek <- krige.control(type.krige = "sk", beta = 5, nugget = 1,
                         cov.model = "exponential", cov.pars = c(1, 0.2))
krige.conv(y.geo, locations = c(0.5, 0.4), krige = simplek)


# sda 104 -----------------------------------------------------------------

library(geoR)
data(ca20)
plot(ca20)

b <- variog(ca20)
plot(b)

par(mfrow=c(1,2))
plot(variog(ca20, estimator.type = "classical"), main = "classical")
plot(variog(ca20, estimator.type = "modulus"), main = "modulus")

b$v

# put 203 - 204------------------------------------------------------------


  # generate data
z <- matrix(rnorm(3000), nrow = 3)
y1 <- 1 + z[1,]
y2 <- 5 + 2 * z[1,] + z[2,]
u <- 2 * (y1 - 1) + z[3,]
miss.index <- which(u < 0)
y2.obs <- y2[-miss.index]
missN <- length(miss.index)

  # imputation bootstrap
y1.use <- y1[miss.index]
B <- 1000 # bootstrap times
sumEY2 <- 0 # init
for (i in 1:B) {
  y1.sample <- sample(y1.use, missN, replace = TRUE)
  y2.boot <- 2 * y1.sample + 3 + rnorm(missN)
  Ey2 <- mean(c(y2.obs, y2.boot))
  sumEY2 = sumEY2 + Ey2
}
sumEY2 / B


library(mice)
impMiss <- mice.impute.norm.boot(y2, 
                                 ry = (u >= 0), 
                                 x = 2 * y1 + 3)
mean(c(y2[-miss.index], impMiss))


N <- 1000
imp <- matrix(nrow = missN, ncol = N)
for (i in 1:N) {
  imp[,i] <- mice.impute.norm(y2,ry = (u >= 0), x = 2 * y1 + 3)  
}

Ey2 <- (colMeans(imp) * missN + sum(y2.obs))/1000

par(mfrow = c(1,2))
plot(colMeans(imp), pch = 20, col = "grey",
     main = "mean of imputation values")
abline(h = mean(colMeans(imp)), col = "red")
plot(Ey2, pch = 20, col = "grey",
     main = "Expectation of Y2")
abline(h = mean(Ey2), col = "red")
abline(h = 5, col = "blue")
