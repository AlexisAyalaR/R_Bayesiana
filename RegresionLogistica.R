#Regresión Logística

library("boot")
data("urine")
?urine
head(urine)

dat = na.omit(urine)

pairs(dat)

library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")


#Selección de variables
X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])
colMeans(X)
apply(X, 2, sd)


#Modelo

#Vamos a usar la distribucion de laplace
ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")


#Jags

library("rjags")

mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)

summary(mod1_sim)

par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))


#Segundo modelo

mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)


#Resultados

dic1
dic2
summary(mod2_sim)

HPDinterval(mod2_csim)

par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))



#Predicciones

(pm_coef = colMeans(mod2_csim))

pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

plot(phat, jitter(dat$r))

(tab0.5 = table(phat > 0.5, data_jags$y))
sum(diag(tab0.5)) / sum(tab0.5)




#Ejercicio
#####
data("PlantGrowth")
?PlantGrowth

boxplot(weight ~ group, data=PlantGrowth)

#JAGS
library("rjags")

#Modelo 1: Misma varianza
#####
mod_string = "model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec)
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0, 1/1e6)
    }
    
    prec ~ dgamma(5/2, 5*1/2)
    sig = sqrt(1/prec)
}"

set.seed(82)
str(PlantGrowth)

data_jags = list(y=PlantGrowth$weight, 
                 grp = as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function(){
  inits = list("mu" = rnorm(3, 0, 100), 
               "prec" = rgamma(1, 1, 1))
}

mod = jags.model(textConnection(mod_string), 
                 data= data_jags, inits = inits, n.chains = 3)
update(mod, 1e3)

mod_sim = coda.samples(model = mod, 
                       variable.names = params,
                       n.iter = 5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim))

plot(mod_sim)

summary(mod_sim)

#Modelo 2: Diferente varianza
#####
mod_string2 = "model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0, 1/1e6)
        prec[j] ~ dgamma(5/2, 5*1/2)
    }
    
    sig = sqrt(1/prec)
}"

set.seed(82)
str(PlantGrowth)

data_jags2 = list(y=PlantGrowth$weight, 
                  grp = as.numeric(PlantGrowth$group))

params2 = c("mu", "sig")

inits2 = function(){
  inits2 = list("mu" = rnorm(3, 0, 100), 
                "prec" = rgamma(3, 1, 1))
}

mod2 = jags.model(textConnection(mod_string2), 
                  data= data_jags2, inits = inits2, n.chains = 3)
update(mod2, 1e3)

mod_sim2 = coda.samples(model = mod2, 
                        variable.names = params2,
                        n.iter = 5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim2))

plot(mod_sim2)

summary(mod_sim2)

#Análisis
#####

dic1 = dic.samples(mod, 1e5)
dic2 = dic.samples(mod2, 1e5)

diff = dic1 - dic2

#Mu3 - Mu1
HPDinterval(mod_csim[,3]-mod_csim[,1])
mean(mod_csim[3]-mod_csim[1])

